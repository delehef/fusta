#![allow(clippy::redundant_field_names)]
use log::*;
use std::ffi::OsStr;
use fuse::*;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::time::{SystemTime, Duration};
use std::fs;
use std::cell::RefCell;
use libc::*;
use maplit::*;

use std::io::SeekFrom;
use std::io::prelude::*;

use fusta::fasta::*;

const TTL: Duration = Duration::from_secs(1);

const FASTA_EXT: &str = "fa";
const SEQ_EXT: &str = "seq";

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
const NO_OVERWRITE: bool = false;

#[derive(Debug, PartialEq, Clone)]
enum FileClass {
    Fasta(std::cell::RefCell<Vec<u8>>),
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

trait VirtualFile {
    fn name(&self) -> &str;
    fn attrs(&self) -> &FileAttr;
    fn mut_attrs(&mut self) -> &mut FileAttr;
    fn ino(&self) -> u64;
    fn class(&self) -> &FileClass;
}


#[derive(Debug)]
struct FragmentFile {
    name: String,
    ino: u64,
    attrs: FileAttr,
    class: FileClass,
}
impl VirtualFile for FragmentFile {
    fn name(&self) -> &str { &self.name }
    fn attrs(&self) -> &FileAttr { &self.attrs }
    fn mut_attrs(&mut self) -> &mut FileAttr { &mut self.attrs }
    fn ino(&self) -> u64 { self.ino }
    fn class(&self) -> &FileClass { &self.class }
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
    fasta_file: FragmentFile,
    seq_file: FragmentFile,
}
impl Fragment {
    fn make_label(id: &str, name: &Option<String>) -> String {
        format!(">{}{}\n", id, name.as_ref().map(|n| format!(" {}", n)).unwrap_or_default())
    }

    fn make_virtual_file(
        ino: u64, name: &str, permissions: u16, size: usize, class: FileClass,
        accessed: SystemTime, modified: SystemTime
    ) -> FragmentFile {
        FragmentFile {
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
                0o664, label.as_bytes().len() + data_size, FileClass::Fasta(RefCell::new(Vec::new())),
                accessed, modified),
            seq_file: Fragment::make_virtual_file(
                seq_ino,
                &format!("{}.{}", id, SEQ_EXT),
                0o664, data_size, FileClass::Seq,
                accessed, modified),
        }
    }

    fn rename(&mut self, new_id: &str) {
        self.id = new_id.to_string();
        self.refresh_virtual_files();
    }

    fn refresh_virtual_files(&mut self) {
        self.fasta_file.name = format!("{}.{}", self.id, FASTA_EXT);
        self.fasta_file.attrs.size = (self.label_size() + self.data_size()) as u64;

        self.seq_file.name = format!("{}.{}", self.id, SEQ_EXT);
        self.seq_file.attrs.size = self.data_size() as u64;
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
                let mut f = fs::File::open(&filename).unwrap_or_else(|_| panic!("Unable to open `{}`", filename));
                f.seek(SeekFrom::Start(*start as u64)).unwrap_or_else(|_| panic!("Unable to seek in `{}`", filename));
                f.read_exact(&mut buffer).unwrap_or_else(|_| panic!("Unable to read from `{}`", filename));
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
                let mut f = fs::File::open(&filename).unwrap_or_else(|_| panic!("Unable to open `{}`", filename));
                f.seek(SeekFrom::Start(*start as u64 + offset as u64)).unwrap_or_else(|_| panic!("Unable to seek in `{}`", filename));
                f.read_exact(&mut buffer).unwrap_or_else(|_| panic!("Unable to read from `{}`", filename));
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

    fn file_from_filename(&self, name: &str) -> Option<&dyn VirtualFile> {
        if self.fasta_file.name == name {
            Some(&self.fasta_file)
        } else if self.seq_file.name == name {
            Some(&self.seq_file)
        } else {
            None
        }
    }

    fn file_from_ino(&self, ino: u64) -> Option<&dyn VirtualFile> {
        if self.fasta_file.ino == ino {
            Some(&self.fasta_file)
        } else if self.seq_file.ino == ino {
            Some(&self.seq_file)
        } else {
            None
        }
    }

    fn mut_file_from_ino(&mut self, ino: u64) -> Option<&mut dyn VirtualFile> {
        if self.fasta_file.ino == ino {
            Some(&mut self.fasta_file)
        } else if self.seq_file.ino == ino {
            Some(&mut self.seq_file)
        } else {
            None
        }
    }
}

pub struct FustaSettings {
    pub mmap: bool,
}

#[derive(Debug)]
struct PendingAppend {
    data: Vec<u8>,
    attrs: FileAttr,
    fasta_ino: u64,
    seq_ino: u64,
}

pub struct FustaFS {
    fragments: Vec<Fragment>,
    metadata: fs::Metadata,
    dir_attrs: BTreeMap<u64, FileAttr>,
    file_attrs: BTreeMap<u64, FileAttr>,
    filename: String,
    settings: FustaSettings,
    current_ino: u64,

    pending_appends: HashMap<String, PendingAppend>,

    info_file_buffer: String,

    dirty: bool,
}


impl FustaFS {
    pub fn new(settings: FustaSettings, filename: &str) -> FustaFS {
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
            dirty: false,
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
        let fasta_file = fs::File::open(filename).unwrap_or_else(|_| panic!("Unable to open `{}`", filename));
        let fragments = FastaReader::new(fasta_file).collect::<Vec<_>>();
        info!("Done.");
        let mut keys = fragments.iter().map(|f| &f.id).collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        if keys.len() != fragments.len() {panic!("Duplicated keys")}

        let file = fs::File::open(filename).unwrap_or_else(|_| panic!("Unable to open `{}`", filename));
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
        if !self.dirty {
            debug!("CONCRETIZE: nothing to do");
            return;
        }

        warn!("========== CONCRETIZING ========");
        let mut index = 0;
        let mut last_start;
        let tmp_filename = format!("{}#fusta#", &self.filename);
        { // Scope to ensure the tmp file is correctly closed
            trace!("Writing fragments");
            let mut tmp_file = fs::File::create(&tmp_filename).unwrap_or_else(|_| panic!("Unable to create `{}`", tmp_filename));
            for fragment in self.fragments.iter_mut() {
                trace!("Writing {}", fragment.id);
                tmp_file.write_all(fragment.label().as_bytes()).unwrap_or_else(|_| panic!("Unable to write to `{}`", tmp_filename));
                index += fragment.label().as_bytes().len();
                last_start = index;
                tmp_file.write_all(&fragment.data()).unwrap_or_else(|_| panic!("Unable to write to `{}`", tmp_filename));
                index += fragment.data().len();

                fragment.data = Backing::File(self.filename.clone(), last_start, index);
                fragment.refresh_virtual_files();
            }
        }
        trace!("Renaming {} to {}", tmp_filename, &self.filename);
        fs::rename(&tmp_filename, &self.filename)
            .unwrap_or_else(|_| panic!("Unable to rename `{}` to `{}`", &tmp_filename, &self.filename));
        warn!("========== DONE ========");
        self.dirty = false;
    }


    fn fragment_from_id(&self, id: &str) -> Option<&Fragment> {
        self.fragments.iter().find(|f| f.id == id)
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
        use ascii_table::*;
        use num_format::*;

        let mut table = AsciiTable::default();
        table.columns.insert(0, Column::with_header("ID"));
        table.columns.insert(1, Column::with_header("Infos"));
        table.columns.insert(2, Column::with_header("Length (bp)"));

        trace!("Making INFO BUFFER");
        let header = format!("{} - {} sequences", &self.filename, self.fragments.len());
        let fragments_infos = self.fragments
            .iter()
            .map(|f| vec![
                f.id.to_owned(),
                f.name.as_ref().unwrap_or(&"".to_string()).to_owned(),
                f.data_size().to_formatted_string(&Locale::en)
            ])
            .collect::<Vec<_>>();
        self.info_file_buffer = format!("{}\n{}\n{}", header, "=".repeat(header.len()), &table.format(&fragments_infos));
        let new_size = self.info_file_buffer.as_bytes().len() as u64;
        self.file_attrs.entry(INFO_FILE).and_modify(|x| x.size = new_size);
    }

    fn is_fasta_file(&self, ino: u64) -> bool {
        self.fragments.iter().any(|f| f.fasta_file.ino == ino)
    }

    fn is_seq_file(&self, ino: u64) -> bool {
        self.fragments.iter().any(|f| f.seq_file.ino == ino)
    }

    fn is_append_file(&self, ino: u64) -> bool {
        self.pending_appends.iter().any(|p| p.1.attrs.ino == ino)
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
                    reply.entry(&TTL, &file.attrs(), 0);
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
            _ => {
                if let Some(file) = self
                    .fragment_from_ino(ino)
                    .and_then(|f| f.file_from_ino(ino))
                {
                    reply.attr(&TTL, &file.attrs())
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
                    match fragment.file_from_ino(ino).expect("No file linked to this fragment").class() {
                        FileClass::Fasta(_) => {
                            let label_size = fragment.label_size() as i64;
                            if offset > label_size as i64 {
                                reply.data(&fragment.chunk(offset - label_size, size))
                            } else {
                                let end = offset as usize + size as usize;
                                let label = fragment.label();
                                let data_chunk = fragment.chunk(0, std::cmp::min(fragment.data_size() as u32, size));
                                match fragment.mut_file_from_ino(ino).expect("No file linked to this fragment").class() {
                                    FileClass::Fasta(header_buffer) => {
                                        if end > header_buffer.borrow().len() {
                                            header_buffer.replace(
                                                label.as_bytes().iter()
                                                    .chain(data_chunk.iter())
                                                    .cloned()
                                                    .take(end as usize)
                                                    .collect::<Vec<_>>()
                                            );
                                        }
                                        reply.data(&header_buffer.borrow()[offset as usize
                                                                           ..
                                                                           std::cmp::min(end, header_buffer.borrow().len())]);
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
                    0          => (FileType::Directory, "..".to_owned()), // TODO
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
                        let length_before = self.fragments.len();
                        self.fragments.retain(|f| f.fasta_file.name != name && f.seq_file.name != name);
                        let length_after = self.fragments.len();
                        // Only mark as dirty if we effectively removed something
                        if length_after != length_before { self.dirty = true; }
                        reply.ok();
                        self.concretize();
                    } else {
                        warn!("UNLINK: unknown file: `{:?}`", name);
                        reply.error(ENOENT);
                    }
            },
            APPEND_DIR => {
                warn!("UNLINK: unauthorized in append virtual dir");
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
                if self.fragments.iter().any(|f| f.id == basename) && NO_OVERWRITE {
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
        if !self.is_writeable(ino) {
            error!("{} is not writeable", ino);
            reply.error(EACCES);
        } else {
            // We write to an existing fragment
            if self.fragment_from_ino(ino).is_some() {
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
                    fragment.refresh_virtual_files();
                    reply.written(data.len() as u32);
                    self.dirty = true;
                } else {
                    panic!("Something went very wrong...")
                }
            }
            // We write to a pending fragment
            else if let Some((name, pending_fragment)) = self
                .pending_appends
                .iter_mut()
                .find(|(_, p)| p.attrs.ino == ino)
            {
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
                    if self.fragment_from_ino(ino).and_then(|f| f.file_from_ino(ino)).unwrap().class().readonly() {
                        reply.error(EACCES);
                    } else {
                        if let Some(file) = self.mut_fragment_from_ino(ino).and_then(|f| f.mut_file_from_ino(ino)) {
                            if let Some(uid) = uid         {file.mut_attrs().uid = uid}
                            if let Some(gid) = gid         {file.mut_attrs().gid = gid}
                            if let Some(atime) = atime     { file.mut_attrs().atime = atime }
                            if let Some(mtime) = mtime     { file.mut_attrs().mtime = mtime }
                            if let Some(chgtime) = chgtime { file.mut_attrs().mtime = chgtime }
                            if let Some(crtime) = crtime   { file.mut_attrs().crtime = crtime }
                            // macOS only
                            if let Some(flags) = flags     { file.mut_attrs().flags = flags }
                            if let Some(mode) = mode       { file.mut_attrs().perm = mode as u16 }
                        }
                        if let Some(size) = size {
                            let size = size as usize;
                            if size == 0 { // Clear the file, called by the truncate syscall
                                if let Some(fragment) = self.mut_fragment_from_ino(ino) {
                                    fragment.data = Backing::Buffer(Vec::new());
                                    fragment.refresh_virtual_files();
                                }
                            } else if size != self.fragment_from_ino(ino).map(Fragment::data_size).unwrap() { // Redim the file
                                let mut fragment = self.mut_fragment_from_ino(ino).unwrap();
                                if !matches!(fragment.data, Backing::Buffer(_)) {
                                    fragment.data = Backing::Buffer(fragment.data().to_vec());
                                }
                                fragment.extend(size);
                                fragment.refresh_virtual_files();
                                self.dirty = true;
                            }
                            self
                                .mut_fragment_from_ino(ino)
                                .and_then(|f| f.mut_file_from_ino(ino))
                                .unwrap().mut_attrs().size = size as u64;
                        }
                        reply.attr(&TTL, &self.fragment_from_ino(ino).and_then(|f| f.file_from_ino(ino)).unwrap().attrs());
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

    fn rename(&mut self, _req: &Request, parent: u64, name: &OsStr, newparent: u64, newname: &OsStr, reply: ReplyEmpty) {
        match parent {
            ROOT_DIR | APPEND_DIR => {
                warn!("RENAME: forbidden here");
                reply.error(EACCES);
            }
            SEQ_DIR | FASTA_DIR => {
                if newparent != parent {
                    error!("Cannot move files out of folder, please copy them instead");
                    reply.error(EACCES);
                } else {
                    let mut new_id = newname.to_str().unwrap().to_string();
                    let newname_path = std::path::Path::new(&new_id);
                    // Remove the artificial extension if it exists
                    if newname_path.extension().map(|ext| ext == SEQ_EXT).unwrap_or(false) {
                        debug!("Using {:?} for {}", newname_path.file_stem(), new_id);
                        new_id = newname_path.file_stem().unwrap().to_str().unwrap().to_string();
                    }
                    // Shortcut if we cannot overwrite existing fragments
                    let replaced_fragment = self.fragment_from_id(&new_id);
                    if replaced_fragment.is_some() && NO_OVERWRITE {
                        error!("Cannot rename {:?} to {}: already existing.", name, new_id);
                        reply.error(EACCES);
                    } else {
                        if replaced_fragment.is_some() {
                            warn!("Replacing {}", new_id);
                            self.fragments.retain(|f| f.id != new_id);
                        }
                        if let Some(ref mut fragment) = self.mut_fragment_from_filename(name.to_str().unwrap()) {
                            fragment.rename(&new_id);
                            info!("Renaming {:?} -> {:?}", name, newname);
                            self.dirty = true;
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
        debug!("RELEASE {}", ino);
        if self.is_writeable(ino) {
            for pending in self.pending_appends.iter() {
                // FS is dirty at the first pending fragment
                self.dirty = true;

                if pending.1.attrs.ino == ino {
                    trace!("RELEASE: {}", pending.0);
                    info!("Dumping...");
                    // Dump the pending in its own file and parse it
                    let mut tmpfile = tempfile::tempfile().expect("Unable to create a temporary file");
                    tmpfile.write_all(&pending.1.data).expect("Unable to write to temporary file");
                    info!("Parsing...");
                    tmpfile.seek(SeekFrom::Start(0)).expect("Unable to seek in temporary file");
                    let fastas = FastaReader::new(&tmpfile).collect::<Vec<_>>();
                    let new_keys = fastas.iter().map(|f| &f.id).collect::<Vec<_>>();
                    let old_keys = self.fragments
                        .iter()
                        .map(|f| &f.id)
                        .cloned()
                        .collect::<Vec<_>>();

                    // Remove the fragments sharing an existing key if overwrite is allowed
                    if !NO_OVERWRITE { self.fragments.retain(|f| !new_keys.contains(&&f.id)) }

                    self.fragments.extend(
                        fastas.iter()
                            .filter_map(|fasta| {
                                if old_keys.contains(&fasta.id) && NO_OVERWRITE {
                                    error!("Skipping `{}`, already existing", &fasta.id);
                                    None
                                } else {
                                    let mut seq = vec![0u8; fasta.pos.1 - fasta.pos.0];
                                    tmpfile.seek(SeekFrom::Start(fasta.pos.0 as u64)).expect("Unable to seek in temporary file");
                                    tmpfile.read_exact(&mut seq).expect("Unable to read from temporary file");

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
        } else {
            debug!("Not a writeable file; ignoring")
        }
        reply.ok();
    }
}
