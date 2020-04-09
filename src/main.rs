#![allow(clippy::redundant_field_names)]

use anyhow::{Context, Result, anyhow, bail};
use fuse::*;
use std::collections::BTreeMap;
use std::collections::HashMap;

use std::ffi::OsStr;
use std::time::{SystemTime, Duration};
use std::sync::mpsc::channel;
use std::fs;
use libc::*;
use daemonize::*;
use maplit::*;

use std::io::SeekFrom;
use std::io::prelude::*;

use clap::*;

use fusta::fasta::*;

use log::*;
use simplelog::*;

const TTL: Duration = Duration::from_secs(1);

const FASTA_EXT: &str = "fa";
const SEQ_EXT: &str = "seq";
const EXTENSIONS: &[&str] = &["FA", "FASTA"];

// Virtual directories
const ROOT_DIR: u64   = 1;
const FASTA_DIR: u64  = 2;
const SEQ_DIR: u64    = 3;
const APPEND_DIR: u64 = 4;

// Pure virtual files
const INFO_FILE: u64  = 10;
const INFO_FILE_NAME: &str = "infos.txt";

// First free ino
const FIRST_INO: u64  = 20;

// Set to true to allow for sequence overwriting
const NO_REPLACE: bool = false;

#[derive(Debug, PartialEq, Clone)]
enum FileClass {
    Fasta(Vec<u8>),
    Seq
}
impl FileClass {
    fn readonly(&self) -> bool {
        match self {
            FileClass::Seq => { false }
            _   => { true }
        }
    }
}

#[derive(Debug)]
struct VirtualFile {
    name: String,
    ino: u64,
    attrs: FileAttr,
    class: FileClass,
}

#[derive(Debug)]
enum Backing {
    File(String, usize, usize),
    Buffer(Vec<u8>),
    MMap(memmap::Mmap)
}
impl Backing {
    fn len(&self) -> usize {
        match self {
            Backing::File(_, start, end) => end - start,
            Backing::Buffer(ref b) => b.len(),
            Backing::MMap(ref mmap) => mmap.len(),
        }
    }
}


#[derive(Debug)]
struct Fragment {
    id: String,
    name: Option<String>,
    data: Backing,
    fasta_file: VirtualFile,
    seq_file: VirtualFile,
}
impl Fragment {
    fn make_label(id: &str, name: &Option<String>) -> String {
        format!(">{}{}\n", id, name.as_ref().map(|n| format!(" {}", n)).unwrap_or_default())
    }

    fn make_virtual_file(
        ino: u64, name: &str, permissions: u16, size: usize, class: FileClass,
        accessed: SystemTime, modified: SystemTime
    ) -> VirtualFile {
        VirtualFile {
            name: name.to_string(),
            ino: ino,
            class: class,
            attrs: FileAttr {
                ino: ino, size: size as u64, blocks: 1,
                atime:  accessed,
                mtime:  modified,
                ctime:  modified,
                crtime: modified,
                kind: FileType::RegularFile, perm: permissions,
                nlink: 0,
                uid: unsafe { libc::geteuid() }, gid: unsafe { libc::getgid() },
                rdev: 0, flags: 0,
            }
        }
    }

    fn new(
        id: &str, name: &Option<String>, data: Backing,
        fasta_ino: u64, seq_ino: u64,
        accessed: SystemTime, modified: SystemTime
    ) -> Fragment {
        let label = Fragment::make_label(id, name);
        let data_size = data.len();
        Fragment {
            id: id.to_string(),
            name: name.clone(),
            data: data,
            fasta_file: Fragment::make_virtual_file(
                fasta_ino,
                &format!("{}.{}", id, FASTA_EXT),
                0o664, label.as_bytes().len() + data_size, FileClass::Fasta(Vec::new()),
                accessed, modified),
            seq_file: Fragment::make_virtual_file(
                seq_ino,
                &format!("{}.{}", id, SEQ_EXT),
                0o664, data_size, FileClass::Seq,
                accessed, modified),
        }
    }

    fn label_size(&self) -> usize {
        self.label().len()
    }

    fn data_size(&self) -> usize {
        match &self.data {
            Backing::File(_, start, end) => end - start,
            Backing::Buffer(ref b)       => b.len(),
            Backing::MMap(ref mmap)      => mmap.len()
        }
    }

    fn total_size(&self) -> usize {
        self.label_size() + self.data_size()
    }

    fn label(&self) -> String {
        Fragment::make_label(&self.id, &self.name)
    }

    fn data(&self) -> Box<[u8]> {
        match &self.data {
            Backing::File(filename, start, _) => {
                let mut buffer = vec![0u8; self.data_size() as usize];
                let mut f = fs::File::open(&filename).expect(&format!("Unable to open `{}`", filename));
                f.seek(SeekFrom::Start(*start as u64)).expect(&format!("Unable to seek in `{}`", filename));
                f.read(&mut buffer).expect(&format!("Unable to read from `{}`", filename));
                buffer.into_boxed_slice()
            },
            Backing::Buffer(ref b) => {
                b[..].into()
            },
            Backing::MMap(ref mmap) => {
                mmap[..].into()
            }
        }
    }

    fn chunk(&self, offset: i64, size: u32) -> Box<[u8]> {
        match &self.data {
            Backing::File(filename, start, _) => {
                let mut buffer = vec![0u8; size as usize];
                let mut f = fs::File::open(&filename).expect(&format!("Unable to open `{}`", filename));
                f.seek(SeekFrom::Start(*start as u64 + offset as u64)).expect(&format!("Unable to seek in `{}`", filename));
                f.read(&mut buffer).expect(&format!("Unable to read from `{}`", filename));
                buffer.into_boxed_slice()
            },
            Backing::Buffer(ref b) => {
                let start = offset as usize;
                let end = std::cmp::min(b.len() as i64, offset + size as i64) as usize;
                b[start .. end].into()
            },
            Backing::MMap(ref mmap) => {
                let start = offset as usize;
                let end = std::cmp::min(mmap.len() as i64, offset + size as i64) as usize;
                mmap[start .. end].into()
            }
        }
    }

    fn extend(&mut self, size: usize) {
        match &mut self.data {
            Backing::Buffer(ref mut b) => {
                b.resize_with(size, Default::default)
            }
            _ => unreachable!()
        }
    }

    fn file_from_filename(&self, name: &str) -> Option<&VirtualFile> {
        if self.fasta_file.name == name {
            Some(&self.fasta_file)
        } else if self.seq_file.name == name {
            Some(&self.seq_file)
        } else {
            None
        }
    }

    fn file_from_ino(&self, ino: u64) -> Option<&VirtualFile> {
        if self.fasta_file.ino == ino {
            Some(&self.fasta_file)
        } else if self.seq_file.ino == ino {
            Some(&self.seq_file)
        } else {
            None
        }
    }

    fn mut_file_from_ino(&mut self, ino: u64) -> Option<&mut VirtualFile> {
        if self.fasta_file.ino == ino {
            Some(&mut self.fasta_file)
        } else if self.seq_file.ino == ino {
            Some(&mut self.seq_file)
        } else {
            None
        }
    }
}

struct FustaSettings {
    mmap: bool,
}

#[derive(Debug)]
struct PendingAppend {
    data: Vec<u8>,
    attrs: FileAttr,
    fasta_ino: u64,
    seq_ino: u64,
}

struct FustaFS {
    fragments: Vec<Fragment>,
    metadata: fs::Metadata,
    dir_attrs: BTreeMap<u64, FileAttr>,
    file_attrs: BTreeMap<u64, FileAttr>,
    filename: String,
    settings: FustaSettings,
    current_ino: u64,

    pending_appends: HashMap<String, PendingAppend>,

    info_file_buffer: String,
}


impl FustaFS {
    fn new(settings: FustaSettings, filename: &str) -> FustaFS {
        let metadata = fs::metadata(filename).unwrap();
        let mut r = FustaFS {
            fragments: Vec::new(),
            filename: String::new(),
            dir_attrs: btreemap! {
                // Virtual folders
                ROOT_DIR   => FustaFS::make_dir_attrs(ROOT_DIR, 0o775),
                SEQ_DIR    => FustaFS::make_dir_attrs(SEQ_DIR, 0o775),
                FASTA_DIR  => FustaFS::make_dir_attrs(FASTA_DIR, 0o555),
                APPEND_DIR => FustaFS::make_dir_attrs(APPEND_DIR, 0o775),
            },
            file_attrs: btreemap! {
                // Pure virtual files
                INFO_FILE  => FustaFS::make_file_attrs(INFO_FILE, 0o444),
            },
            metadata: metadata,
            settings: settings,
            current_ino: FIRST_INO,
            pending_appends: HashMap::new(),
            info_file_buffer: String::new(),
        };

        r.read_fasta(filename);
        r
    }

    fn make_dir_attrs(ino: u64, perms: u16) -> FileAttr {
        FileAttr {
            ino: ino,
            size: 0,
            blocks: 0,
            atime: std::time::SystemTime::now(),
            mtime: std::time::SystemTime::now(),
            ctime: std::time::SystemTime::now(),
            crtime: std::time::SystemTime::now(),
            kind: FileType::Directory,
            perm: perms,
            nlink: 1,
            uid: unsafe { libc::geteuid() },
            gid: unsafe { libc::getgid() },
            rdev: 0,
            flags: 0,
        }
    }

    fn make_file_attrs(ino: u64, perms: u16) -> FileAttr {
        FileAttr {
            ino: ino,
            size: 0,
            blocks: 0,
            atime: std::time::SystemTime::now(),
            mtime: std::time::SystemTime::now(),
            ctime: std::time::SystemTime::now(),
            crtime: std::time::SystemTime::now(),
            kind: FileType::RegularFile,
            perm: perms,
            nlink: 1,
            uid: unsafe { libc::geteuid() },
            gid: unsafe { libc::getgid() },
            rdev: 0,
            flags: 0,
        }
    }

    fn new_ino(&mut self) -> u64 {
        let r = self.current_ino;
        self.current_ino += 1;
        trace!("New ino: {}", r);
        r
    }

    fn read_fasta(&mut self, filename: &str) {
        info!("Reading {}...", filename);
        let fasta_file = fs::File::open(filename).expect(&format!("Unable to open `{}`", filename));
        let fragments = FastaReader::new(fasta_file).collect::<Vec<_>>();
        info!("Done.");
        let mut keys = fragments.iter().map(|f| &f.id).collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        if keys.len() != fragments.len() {panic!("Duplicated keys")}

        let file = fs::File::open(filename).expect(&format!("Unable to open `{}`", filename));
        self.filename = filename.to_owned();

        self.fragments = fragments.iter()
            .map(|fragment| {
                Fragment::new(
                    &fragment.id, &fragment.name,
                    if self.settings.mmap {
                        Backing::MMap(unsafe {
                            memmap::MmapOptions::new()
                                .offset(fragment.pos.0 as u64)
                                .len(fragment.len)
                                .map(&file)
                                .unwrap()
                        })
                    } else {
                        Backing::File(filename.to_owned(), fragment.pos.0, fragment.pos.1)
                    },
                    self.new_ino(), self.new_ino(),
                    self.metadata.accessed().unwrap(),
                    self.metadata.modified().unwrap(),
                )
            })
            .collect::<Vec<_>>();
        self.make_info_buffer();
    }

    fn concretize(&mut self) {
        let mut index = 0;
        let mut last_start;
        let tmp_filename = format!("{}#fusta#", &self.filename);
        { // Scope to ensure the tmp file is correctly closed
            trace!("Writing fragments");
            let mut tmp_file = fs::File::create(&tmp_filename).expect(&format!("Unable to create `{}`", tmp_filename));
            for fragment in self.fragments.iter_mut() {
                trace!("Writing {}", fragment.id);
                tmp_file.write_all(fragment.label().as_bytes()).expect(&format!("Unable to write to `{}`", tmp_filename));
                index += fragment.label().as_bytes().len();
                last_start = index;
                tmp_file.write_all(&fragment.data()).expect(&format!("Unable to write to `{}`", tmp_filename));
                index += fragment.data().len();

                fragment.data = Backing::File(self.filename.clone(), last_start, index);
            }
        }
        trace!("Renaming {} to {}", tmp_filename, &self.filename);
        fs::rename(&tmp_filename, &self.filename).expect(&format!("Unable to rename `{}` to `{}`", &tmp_filename, &self.filename));
        error!("========== CONCRETIZING ========")
    }

    fn fragment_from_ino(&self, ino: u64) -> Option<&Fragment> {
        self.fragments.iter().find(|f| f.fasta_file.ino == ino || f.seq_file.ino == ino)
    }

    fn mut_fragment_from_ino(&mut self, ino: u64) -> Option<&mut Fragment> {
        self.fragments.iter_mut().find(|f| f.fasta_file.ino == ino || f.seq_file.ino == ino)
    }

    fn fragment_from_filename(&self, name: &str) -> Option<&Fragment> {
        self.fragments.iter().find(|f| f.fasta_file.name == name || f.seq_file.name == name)
    }

    fn mut_fragment_from_filename(&mut self, name: &str) -> Option<&mut Fragment> {
        self.fragments.iter_mut().find(|f| f.fasta_file.name == name || f.seq_file.name == name)
    }

    fn make_info_buffer(&mut self) {
        error!("Making INFO BUFFER");
        let header = format!("{} - {} sequences", &self.filename, self.fragments.len());
        let fragments_infos = self.fragments
            .iter()
            .map(|f| format!("{} {} - {} bp", f.id, f.name.as_ref().unwrap_or(&"".to_string()), f.data_size()))
            .collect::<Vec<_>>();
        self.info_file_buffer = format!("{}\n{}\n{}\n", header, "=".repeat(header.len()), &fragments_infos.join("\n"));
        self.file_attrs.get_mut(&INFO_FILE).unwrap().size = self.info_file_buffer.as_bytes().len() as u64;
    }

    fn is_fasta_file(&self, ino: u64) -> bool {
        self.fragments.iter().find(|f| f.fasta_file.ino == ino).is_some()
    }

    fn is_seq_file(&self, ino: u64) -> bool {
        self.fragments.iter().find(|f| f.seq_file.ino == ino).is_some()
    }

    fn is_append_file(&self, ino: u64) -> bool {
        self.pending_appends.iter().find(|p| p.1.attrs.ino == ino).is_some()
    }

    fn is_writeable(&self, ino: u64) -> bool {
        self.is_append_file(ino) || self.is_seq_file(ino)
    }
}


impl Filesystem for FustaFS {
    fn lookup(&mut self, _req: &Request, parent: u64, name: &OsStr, reply: ReplyEntry) {
        let name = name.to_str().unwrap();
        match parent {
            ROOT_DIR => {
                match name {
                    "fasta" => {
                        reply.entry(&TTL, &self.dir_attrs[&FASTA_DIR], 0);
                    },
                    "seqs" => {
                        reply.entry(&TTL, &self.dir_attrs[&SEQ_DIR], 0);
                    },
                    "append" => {
                        reply.entry(&TTL, &self.dir_attrs[&APPEND_DIR], 0);
                    },
                    INFO_FILE_NAME => {
                        reply.entry(&TTL, &self.file_attrs[&INFO_FILE], 0);
                    }
                    _ => {
                        reply.error(ENOENT);
                    }
                }
            },
            SEQ_DIR | FASTA_DIR => {
                if let Some(file) = self
                    .fragment_from_filename(name)
                    .and_then(|f| f.file_from_filename(name))
                {
                    reply.entry(&TTL, &file.attrs, 0);
                } else {
                    reply.error(ENOENT);
                }
            },
            _ => {
                warn!("LOOKUP: parent {} does not exist", parent);
                reply.error(ENOENT);
            }
        }
    }

    fn getattr(&mut self, _req: &Request, ino: u64, reply: ReplyAttr) {
        match ino {
            ROOT_DIR | SEQ_DIR | APPEND_DIR | FASTA_DIR => {
                reply.attr(&TTL, &self.dir_attrs.get(&ino).unwrap())
            }
            INFO_FILE => {
                reply.attr(&TTL, &self.file_attrs.get(&ino).unwrap())
            }
            _          => {
                if let Some(file) = self
                    .fragment_from_ino(ino)
                    .and_then(|f| f.file_from_ino(ino))
                {
                    reply.attr(&TTL, &file.attrs)
                } else {
                    warn!("GETATTR: ino `{}` does not exist", ino);
                    reply.error(ENOENT)
                }
            },
        }
    }

    fn read(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, size: u32, reply: ReplyData) {
        debug!("READING {}", ino);
        match ino {
            INFO_FILE => {
                self.make_info_buffer();
                let start = offset as usize;
                let end = std::cmp::min(start + size as usize, self.info_file_buffer.len());
                reply.data(&self.info_file_buffer.as_bytes()[start .. end]);
            },
            _ => {
                if self.fragment_from_ino(ino).is_some() {
                    let fragment = match self.mut_fragment_from_ino(ino) {
                        Some(f) => f,
                        _ => {
                            error!("No fragment linked to ino {}", ino);
                            return;
                        }
                    };
                    match fragment.file_from_ino(ino).expect("No file linked to this fragment").class {
                        FileClass::Fasta(_) => {
                            let label_size = fragment.label_size() as i64;
                            if offset > label_size as i64 {
                                reply.data(&fragment.chunk(offset - label_size, size))
                            } else {
                                let end = offset as usize + size as usize;
                                let label = fragment.label();
                                let data_chunk = fragment.chunk(0, std::cmp::min(fragment.data_size() as u32, size));
                                match fragment.mut_file_from_ino(ino).expect("No file linked to this fragment").class {
                                    FileClass::Fasta(ref mut header_buffer) => {
                                        if end > header_buffer.len() {
                                            *header_buffer = label.as_bytes().iter()
                                                .chain(data_chunk.iter())
                                                .cloned()
                                                .take(end as usize)
                                                .collect::<Vec<_>>();
                                        }
                                        reply.data(&header_buffer[offset as usize .. std::cmp::min(end, header_buffer.len())]);
                                    }
                                    _ => { panic!("WTF") }
                                }
                            }
                        }
                        FileClass::Seq => {
                            reply.data(&fragment.chunk(offset, size))
                        }
                    }
                } else {
                    warn!("READ: {} is not a file", ino);
                    reply.error(ENOENT);
                }
            }
        }
    }

    fn readdir(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, mut reply: ReplyDirectory) {
        match ino {
            ROOT_DIR => {
                let entries = btreemap! {
                    ROOT_DIR   => (FileType::Directory, ".".to_owned()),
                    45         => (FileType::Directory, "..".to_owned()), // TODO
                    FASTA_DIR  => (FileType::Directory, "fasta".to_owned()),
                    SEQ_DIR    => (FileType::Directory, "seqs".to_owned()),
                    APPEND_DIR => (FileType::Directory, "append".to_owned()),
                    INFO_FILE  => (FileType::RegularFile, INFO_FILE_NAME.to_owned()),
                };
                for (o, (ino, entry)) in entries.iter().enumerate().skip(offset as usize) {
                    reply.add(*ino, o as i64 + 1, entry.0, entry.1.to_owned());
                }
                reply.ok();
            },
            FASTA_DIR => {
                let mut entries = btreemap! {
                    FASTA_DIR  => (FileType::Directory, ".".to_owned()),
                    ROOT_DIR => (FileType::Directory, "..".to_owned()),
                };
                for f in self.fragments.iter() {
                    entries.insert(f.fasta_file.ino, (FileType::RegularFile, f.fasta_file.name.clone()));
                }

                for (o, (ino, entry)) in entries.iter().enumerate().skip(offset as usize) {
                    reply.add(*ino, o as i64 + 1, entry.0, entry.1.to_owned());
                }
                reply.ok();
            },
            SEQ_DIR => {
                let mut entries = btreemap! {
                    SEQ_DIR  => (FileType::Directory, ".".to_owned()),
                    ROOT_DIR => (FileType::Directory, "..".to_owned()),
                };
                for f in self.fragments.iter() {
                    entries.insert(f.seq_file.ino, (FileType::RegularFile, f.seq_file.name.clone()));
                }

                for (o, (ino, entry)) in entries.iter().enumerate().skip(offset as usize) {
                    reply.add(*ino, o as i64 + 1, entry.0, entry.1.to_owned());
                }
                reply.ok();
            },
            APPEND_DIR => {
                let entries = btreemap! {
                    APPEND_DIR => (FileType::Directory, ".".to_owned()),
                    ROOT_DIR   => (FileType::Directory, "..".to_owned()),
                };
                for (o, (ino, entry)) in entries.iter().enumerate().skip(offset as usize) {
                    reply.add(*ino, o as i64 + 1, entry.0, entry.1.to_owned());
                }
                reply.ok();
            }
            _ => {
                warn!("{} is not a directory", ino);
                reply.error(ENOENT);
            }
        }

    }

    fn unlink(&mut self, _req: &Request, parent: u64, name: &OsStr, reply: ReplyEmpty) {
        match parent {
            ROOT_DIR => {
                warn!("UNLINK: cannot remove `{:?}`", name);
            },
            SEQ_DIR | FASTA_DIR => {
                let name = name.to_str().unwrap();
                if self
                    .fragment_from_filename(name)
                    .and_then(|f| f.file_from_filename(name))
                    .is_some() {
                        self.fragments.retain(|f| f.fasta_file.name != name && f.seq_file.name != name);
                        reply.ok();
                        self.concretize();
                    } else {
                        warn!("UNLINK: unknown file: `{:?}`", name);
                        reply.error(ENOENT);
                    }
            },
            APPEND_DIR => {
                warn!("UNLINK: cannot remove from append virtual dir");
                reply.error(EACCES);
            }
            _ => {
                warn!("UNLINK: parent {} does not exist", parent);
                reply.error(ENOENT);
            }
        }
    }

    fn mknod(&mut self, _req: &Request, parent: u64, name: &OsStr, _mode: u32, _rdev: u32, reply: ReplyEntry) {
        match parent {
            ROOT_DIR | SEQ_DIR | FASTA_DIR => {
                warn!("MKNOD: writing in {} is forbidden", parent);
                reply.error(EACCES);
            },
            APPEND_DIR => {
                let name = name.to_str().unwrap();
                let basename = std::path::Path::new(name).file_stem().unwrap().to_str().unwrap();
                // If pathname already exists [...], this call fails with an EEXIST error.
                if self.fragments.iter().any(|f| f.id == basename) && NO_REPLACE {
                    error!("Cannot create `{:?}`, already exists", name);
                    reply.error(EEXIST);
                    return
                }

                trace!("Creating pending append {}", name);
                let attrs = FileAttr {
                    ino: self.new_ino(),
                    size: 0,
                    blocks: 0,
                    atime: std::time::SystemTime::now(),
                    mtime: std::time::SystemTime::now(),
                    ctime: std::time::SystemTime::now(),
                    crtime: std::time::SystemTime::now(),
                    kind: FileType::RegularFile,
                    perm: 0o775,
                    nlink: 0,
                    uid: unsafe { libc::geteuid() },
                    gid: unsafe { libc::getgid() },
                    rdev: 0,
                    flags: 0,
                };
                let pending = PendingAppend {
                    data: Vec::new(),
                    attrs: attrs,
                    seq_ino: self.new_ino(),
                    fasta_ino: self.new_ino(),
                };
                reply.entry(&TTL, &pending.attrs, 0);
                self.pending_appends.insert(basename.to_string(), pending);
            },
            _ => {
                warn!("MKNOD: parent {} does not exist", parent);
                reply.error(ENOENT);
            }
        }
    }

    fn write(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, data: &[u8], _flags: u32, reply: ReplyWrite) {
        if self.fragment_from_ino(ino).is_some() { // We write to an existing fragment
            if self.fragment_from_ino(ino).and_then(|f| f.file_from_ino(ino)).unwrap().class.readonly() {
                reply.error(EACCES);
            } else {
                // If it's a raw seq file, we can write
                let mut fragment = self.mut_fragment_from_ino(ino).expect("Something went very wrong");

                // As soon as there's a write, we have to switch to a buffer-backed storage
                if !matches!(fragment.data, Backing::Buffer(_)) {
                    fragment.data = Backing::Buffer(fragment.data().to_vec());
                }

                // Ensure that the backing buffer is big enough
                let max_size = offset as usize + data.len();
                if max_size > fragment.data_size() { fragment.extend(max_size) }

                // Then finally write the data
                if let Backing::Buffer(b) = &mut fragment.data {
                    let start = offset as usize;
                    let end = start + data.len();
                    b.splice(start .. end, data.iter().cloned());
                    reply.written(data.len() as u32);
                } else {
                    panic!("Something went very wrong...")
                }
            }
        } else if let Some((name, pending_fragment)) = self
            .pending_appends
            .iter_mut()
            .find(|(_, p)| p.attrs.ino == ino)
        { // We write to a pending fragment
            let start = offset as usize;
            let end = start + data.len();
            if end > pending_fragment.data.len() {
                pending_fragment.data.resize_with(end, Default::default);
            }
            pending_fragment.data.splice(start .. end, data.iter().cloned());
            reply.written(data.len() as u32);
            trace!("\tWriting to {}", name);
        } else {
            reply.error(ENOENT);
        }
    }

    fn setattr(&mut self, _req: &Request<'_>, ino: u64,
               mode: Option<u32>, uid: Option<u32>, gid: Option<u32>,
               size: Option<u64>,
               atime: Option<SystemTime>, mtime: Option<SystemTime>,
               _fh: Option<u64>,
               crtime: Option<SystemTime>, chgtime: Option<SystemTime>, bkuptime: Option<SystemTime>,
               flags: Option<u32>,
               reply: ReplyAttr) {
        trace!("SETATTR");
        trace!("mode       {:?}", mode);
        trace!("gid        {:?}", gid);
        trace!("uid        {:?}", uid);
        trace!("size       {:?}", size);
        trace!("atime      {:?}", atime);
        trace!("mtime      {:?}", mtime);
        trace!("bkuptime   {:?}", bkuptime);
        trace!("chgtime    {:?}", chgtime);
        trace!("crtime     {:?}", crtime);
        trace!("flags      {:?}", flags);

        match ino {
            ROOT_DIR | SEQ_DIR | FASTA_DIR => { reply.error(EACCES) }
            INFO_FILE => { reply.error(EACCES) }
            _ => {
                if self.fragment_from_ino(ino).is_some() {
                    if self.fragment_from_ino(ino).and_then(|f| f.file_from_ino(ino)).unwrap().class.readonly() {
                        reply.error(EACCES);
                    } else {
                        if let Some(file) = self.mut_fragment_from_ino(ino).and_then(|f| f.mut_file_from_ino(ino)) {
                            if let Some(uid) = uid         {file.attrs.uid = uid}
                            if let Some(gid) = gid         {file.attrs.gid = gid}
                            if let Some(atime) = atime     { file.attrs.atime = atime }
                            if let Some(mtime) = mtime     { file.attrs.mtime = mtime }
                            if let Some(chgtime) = chgtime { file.attrs.mtime = chgtime }
                            if let Some(crtime) = crtime   { file.attrs.crtime = crtime }
                            // macOS only
                            if let Some(flags) = flags     { file.attrs.flags = flags }
                            if let Some(mode) = mode       { file.attrs.perm = mode as u16 }
                        }
                        if let Some(size) = size {
                            let size = size as usize;
                            if size == 0 { // Clear the file, called by the truncate syscall
                                if let Some(f) = self.mut_fragment_from_ino(ino) {
                                    f.data = Backing::Buffer(Vec::new());
                                }
                            } else if size != self.fragment_from_ino(ino).map(Fragment::data_size).unwrap() { // Redim the file
                                let mut fragment = self.mut_fragment_from_ino(ino).unwrap();
                                if !matches!(fragment.data, Backing::Buffer(_)) {
                                    fragment.data = Backing::Buffer(fragment.data().to_vec());
                                }
                                fragment.extend(size)
                            }
                            self
                                .mut_fragment_from_ino(ino)
                                .and_then(|f| f.mut_file_from_ino(ino))
                                .unwrap().attrs.size = size as u64;
                        }
                        reply.attr(&TTL, &self.fragment_from_ino(ino).and_then(|f| f.file_from_ino(ino)).unwrap().attrs);
                    }
                } else if let Some((name, pending_fragment)) = self
                    .pending_appends
                    .iter_mut()
                    .find(|(_, p)| p.attrs.ino == ino)
                {
                    if let Some(size) = size {
                        trace!("\tResizing {} @{}", name, size);
                        pending_fragment.data.resize_with(size as usize, Default::default);
                        reply.attr(&TTL, &pending_fragment.attrs);
                    }
                } else {
                    warn!("\t{:?} does not exist", ino);
                    reply.error(ENOENT);
                }
            }
        }
    }


    fn destroy(&mut self, _req: &Request) {
        self.concretize()
    }

    /// Rename a file.
    fn rename(&mut self, _req: &Request, parent: u64, name: &OsStr, newparent: u64, newname: &OsStr, reply: ReplyEmpty) {
        match parent {
            ROOT_DIR | APPEND_DIR => {
                warn!("RENAME: forbidden here");
                reply.error(EACCES);
            }
            SEQ_DIR | FASTA_DIR => {
                if newparent != parent {
                    error!("RENAME: Cannot move files oout of folder, please copy them");
                    reply.error(EACCES);
                } else {
                    let keys = self.fragments.iter().map(|f| f.id.clone()).collect::<Vec<_>>();
                    let mut newname = newname.to_str().unwrap().to_string();
                    let newname_path = std::path::Path::new(&newname);
                    // Remove the artificial extension
                    if newname_path.extension().map(|ext| ext == SEQ_EXT).unwrap_or(false) {
                        trace!("Using {:?} for {}", newname_path.file_stem(), newname);
                        newname = newname_path.file_stem().unwrap().to_str().unwrap().to_string();
                    }
                    if keys.contains(&newname) && NO_REPLACE {
                        error!("Cannot rename {:?} to {}: already existing.", name, newname);
                        reply.error(EACCES);
                    } else {
                        if keys.contains(&newname) {
                            warn!("Replacing {}", newname);
                            self.fragments.retain(|f| f.id != newname);
                        }
                        if let Some(ref mut fragment) = self.mut_fragment_from_filename(name.to_str().unwrap()) {
                            fragment.id = newname.to_string();
                            info!("Renaming {:?} -> {:?}", name, newname);
                            self.concretize();
                            reply.ok()
                        } else {
                            warn!("\t{:?} does not exist", name);
                            reply.error(ENOENT);
                        }
                    }
                }
            }
            _ => {
                warn!("RENAME: unknown parent {}", parent);
                reply.error(ENOENT);
            }
        }
    }

    fn fsync(&mut self, _req: &Request, _ino: u64, _fh: u64, _datasync: bool, reply: ReplyEmpty) {
        trace!("FSYNC");
        self.concretize();
        reply.ok();
    }

    fn fsyncdir (&mut self, _req: &Request, _ino: u64, _fh: u64, _datasync: bool, reply: ReplyEmpty) {
        trace!("FSYNCDIR");
        self.concretize();
        reply.ok();
    }

    fn flush(&mut self, _req: &Request, _ino: u64, _fh: u64, _lock_owner: u64, reply: ReplyEmpty) {
        trace!("FLUSH");
        self.concretize();
        reply.ok();
    }

    fn release(&mut self, _req: &Request, ino: u64, _fh: u64, _flags: u32, _lock_owner: u64, _flush: bool, reply: ReplyEmpty) {
        trace!("RELEASE {}", ino);
        match ino {
            _ if self.is_writeable(ino) => {
                for pending in self.pending_appends.iter() {
                    if pending.1.attrs.ino == ino {
                        trace!("RELEASE: {}", pending.0);
                        info!("Dumping...");
                        let mut tmpfile = tempfile::tempfile().expect(&format!("Unable to create a temporary file"));
                        tmpfile.write_all(&pending.1.data).expect(&format!("Unable to write to temporary file"));

                        info!("Parsing...");
                        tmpfile.seek(SeekFrom::Start(0)).expect(&format!("Unable to seek in temporary file"));
                        let fastas = FastaReader::new(&tmpfile).collect::<Vec<_>>();

                        let new_keys = fastas.iter().map(|f| &f.id).collect::<Vec<_>>();
                        let old_keys = self.fragments
                            .iter()
                            .map(|f| &f.id)
                            .cloned()
                            .collect::<Vec<_>>();

                        if !NO_REPLACE { self.fragments.retain(|f| !new_keys.contains(&&f.id)) }

                        self.fragments.extend(
                            fastas.iter()
                                .filter_map(|fasta| {
                                    let mut seq = vec![0u8; fasta.pos.1 - fasta.pos.0];
                                    tmpfile.seek(SeekFrom::Start(fasta.pos.0 as u64)).expect(&format!("Unable to seek in temporary file"));
                                    let _ = tmpfile.read(&mut seq).unwrap();

                                    if old_keys.contains(&fasta.id) && NO_REPLACE {
                                        error!("Skipping {}, already existing", &fasta.id);
                                        None
                                    } else {
                                        Some(Fragment::new(
                                            &fasta.id, &fasta.name,
                                            Backing::Buffer(seq),
                                            pending.1.seq_ino, pending.1.fasta_ino,
                                            pending.1.attrs.atime, pending.1.attrs.mtime
                                        ))
                                    }
                                }))
                    }
                }
                self.concretize();
            }
            _ => { debug!("Not a writeable file; ignoring") }
        }
        reply.ok();
    }

    /// Get file system statistics.
    fn statfs(&mut self, _req: &Request, _ino: u64, reply: ReplyStatfs) {
        trace!("STATFS");
        reply.statfs(
            0,                           // blocks
            0,                           // bfree
            0,                           // bavail
            self.fragments.len() as u64, // files
            0,                           // ffree
            512,                         //bsize
            255,                         // namelen
            0                            // frsize
        );
    }
}

#[derive(Debug, Clone)]
struct RunEnvironment {
    mountpoint: std::path::PathBuf,
    created_mountpoint: bool,
}
fn main() -> Result<()> {
    let args = App::new("fusta")
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .version(crate_version!())
        .author(crate_authors!())
        .arg(Arg::with_name("FASTA")
             .help("A (multi)FASTA file containing the sequences to mount")
             .required(true)
             .index(1))
        .arg(Arg::with_name("mountpoint")
             .short("o")
             .long("mountpoint")
             .help("Specifies the directory to use as mountpoint")
             .default_value("fusta")
             .takes_value(true)
             .required(true))

        .arg(Arg::with_name("v")
             .short("v")
             .multiple(true)
             .help("Sets the level of verbosity"))
        .arg(Arg::with_name("daemon")
             .short("D")
             .long("daemon")
             .help("Launch in the background; will automatically quit when unmounted"))

        .arg(Arg::with_name("mmap")
             .short("M")
             .long("mmap")
             .help("Use mmap instead of seek to browse FASTA fragments. Faster, but memory hungry"))

        .arg(Arg::with_name("nonempty")
             .short("E")
             .long("non-empty")
             .help("Perform the mount even if the destination folder is not empty"))
        .get_matches();

    let log_level = match args.occurrences_of("v") {
        0 => LevelFilter::Warn,
        1 => LevelFilter::Info,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Trace,
    };
    let log_config = ConfigBuilder::new()
        .set_time_format_str("")
        .build();

    if args.is_present("daemon") {
        let mut log_file = std::env::temp_dir();
        log_file.push("fusta.log");
        WriteLogger::new(log_level, log_config, std::fs::File::create(log_file)?);

        let mut pid_file = std::env::temp_dir();
        pid_file.push("fusta.pid");

        Daemonize::new()
            .pid_file(pid_file)
            .working_directory(std::env::current_dir().expect("Unable to read current directory"))
            .start()?;
    } else {
        TermLogger::init(log_level, log_config, TerminalMode::Mixed).expect("Unable to initialize logger");
    }



    let fasta_file = value_t!(args, "FASTA", String)?;
    let mountpoint = value_t!(args, "mountpoint", String)?;
    let mut fuse_options: Vec<&OsStr> = vec![
        &OsStr::new("-o"), &OsStr::new("auto_unmount"),
        &OsStr::new("-o"), &OsStr::new("default_permissions"),
    ];
    let settings = FustaSettings {
        mmap: args.is_present("mmap"),
    };

    info!("Using MMAP:      {}", settings.mmap);

    let fs = FustaFS::new(settings, &fasta_file);
    let mut env = RunEnvironment {
        mountpoint: std::path::PathBuf::from(mountpoint),
        created_mountpoint: false,
    };

    if !env.mountpoint.exists() {
        std::fs::create_dir(&env.mountpoint)?;
        env.created_mountpoint = true;
    }
    if !env.mountpoint.is_dir() {
        bail!("mount point `{:?}` is not a directory", env.mountpoint);
    }

    if args.is_present("nonempty") {
        fuse_options.push(&OsStr::new("-o"));
        fuse_options.push(&OsStr::new("nonempty"));
    }
    if !args.is_present("nonempty") && std::fs::read_dir(&env.mountpoint)?.take(1).count() != 0 {
        bail!(
            "mount point {:?} is not empty, use the -E flag to mount in a non-empty directory.",
            env.mountpoint
        );
    }

    let (tx, rx) = channel();


    {
        let tx_ctrlc = tx.clone();
        ctrlc::set_handler(move || {
            info!("Ctrl-C received, exiting.");
            let _ = tx_ctrlc.send(());
        })?;
    }

    {
        let tx_fuse = tx.clone();
        let env = env.clone();
        let _ = std::thread::spawn(move || {
            match fuse::mount(fs, &env.mountpoint, &fuse_options) {
                Ok(()) => {},
                _      => { error!("Unable to mount the FUSE filesystem"); std::process::exit(1); }
            }
            let _ = tx_fuse.send(());
        });
    }

    rx.recv()?;
    cleanup(&env)?;

    Ok(())
}

fn cleanup(env: &RunEnvironment) -> Result<()> {
    if env.created_mountpoint {
        warn!("You can now safely remove the {:?} directory", env.mountpoint);
    }
    Ok(())
}
