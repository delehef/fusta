use fuse::*;
use std::collections::BTreeMap;

use std::ffi::OsStr;
use std::time::{SystemTime, Duration};
use std::fs;
use libc::*;
use daemonize::*;

use memmap;
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

const ROOT_DIR: u64 = 1;
const FASTA_DIR: u64 = 2;
const SEQ_DIR: u64 = 3;
const FIRST_INO: u64 = 5;

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

    fn with_empty_buffer(
        id: &str, fasta_ino: u64, seq_ino: u64,
        accessed: SystemTime, modified: SystemTime
    ) -> Fragment {
        let label = format!(">{}\n", id);
        Fragment {
            id: id.to_string(),
            name: None,
            data: Backing::Buffer(Vec::new()),
            fasta_file: Fragment::make_virtual_file(
                fasta_ino,
                &format!("{}.{}", id, FASTA_EXT),
                0o444, label.as_bytes().len(), FileClass::Fasta(Vec::new()),
                accessed, modified),
            seq_file: Fragment::make_virtual_file(
                seq_ino,
                &format!("{}.{}", id, SEQ_EXT),
                0o664, 0, FileClass::Seq,
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
                let mut f = fs::File::open(&filename).unwrap();
                f.seek(SeekFrom::Start(*start as u64)).unwrap();
                let _ = f.read(&mut buffer).unwrap();
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
                let mut f = fs::File::open(&filename).unwrap();
                f.seek(SeekFrom::Start(*start as u64 + offset as u64)).unwrap();
                let _ = f.read(&mut buffer).unwrap();
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

    fn to_filename(&self) -> String {
        format!("{}", self.id)
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

struct FustaFS {
    fragments: Vec<Fragment>,
    metadata: fs::Metadata,
    dir_attrs: BTreeMap<u64, FileAttr>,
    filename: String,
    settings: FustaSettings,
    current_ino: u64,
}


impl FustaFS {
    fn new(settings: FustaSettings, filename: &str) -> FustaFS {
        let metadata = fs::metadata(filename).unwrap();
        let mut r = FustaFS {
            fragments: Vec::new(),
            filename: String::new(),
            dir_attrs: maplit::btreemap! {
                ROOT_DIR => FustaFS::make_dir_attrs(ROOT_DIR, 0o775),
                SEQ_DIR => FustaFS::make_dir_attrs(SEQ_DIR, 0o775),
                FASTA_DIR => FustaFS::make_dir_attrs(FASTA_DIR, 0o555),
            },
            metadata: metadata,
            settings: settings,
            current_ino: FIRST_INO,
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

    fn new_ino(&mut self) -> u64 {
        let r = self.current_ino;
        self.current_ino += 1;
        debug!("New ino: {}", r);
        r
    }

    fn read_fasta(&mut self, filename: &str) {
        info!("Reading {}...", filename);
        let fasta_file = fs::File::open(filename).unwrap();
        let fragments = FastaReader::new(fasta_file, false).collect::<Vec<_>>();
        info!("Done.");
        let mut keys = fragments.iter().map(|f| &f.id).collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        if keys.len() != fragments.len() {panic!("Duplicated keys")}

        let file = fs::File::open(filename).unwrap();
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
    }

    /// if force is set to true, we still force the rewrite so that "phantom" (deleted) fragments
    /// are effectively removed.
    /// Otherwise, the memory-backed fragments are written down on disk.
    fn concretize(&mut self, force: bool) {
        if !force && !self.fragments.iter().any(|f| matches!(f.data, Backing::Buffer(_))) {
            debug!("Nothing to concretize; leaving");
            return;
        }

        let tmp_filename = format!("{}#fusta#", &self.filename);
        { // Scope to ensure the tmp file is correctly closed
            trace!("Writing fragments");
            let mut tmp_file = fs::File::create(&tmp_filename).unwrap();
            for fragment in self.fragments.iter() {
                trace!("Writing {}", fragment.id);
                tmp_file.write_all(fragment.label().as_bytes()).unwrap();
                tmp_file.write_all(&fragment.data()).unwrap();
            }
        }
        trace!("Renaming {} to {}", tmp_filename, &self.filename);
        fs::rename(tmp_filename, &self.filename).unwrap();
        trace!("Rebuilding index");
        let filename = self.filename.clone();
        self.read_fasta(&filename);
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
            ROOT_DIR => reply.attr(&TTL, &self.dir_attrs.get(&ROOT_DIR).unwrap()),
            SEQ_DIR => reply.attr(&TTL, &self.dir_attrs.get(&SEQ_DIR).unwrap()),
            FASTA_DIR => reply.attr(&TTL, &self.dir_attrs.get(&FASTA_DIR).unwrap()),
            _ => {
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
        if self.fragment_from_ino(ino).is_some() {
            let fragment = self.mut_fragment_from_ino(ino).unwrap();
            match fragment.file_from_ino(ino).unwrap().class {
                FileClass::Fasta(_) => {
                    let label_size = fragment.label_size() as i64;
                    if offset > label_size as i64 {
                        reply.data(&fragment.chunk(offset - label_size, size))
                    } else {
                        let end = offset as usize + size as usize;
                        let label = fragment.label().clone();
                        let data_chunk = fragment.chunk(0, std::cmp::min(fragment.data_size() as u32, size));
                        match fragment.mut_file_from_ino(ino).unwrap().class {
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

    fn readdir(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, mut reply: ReplyDirectory) {
        match ino {
            ROOT_DIR => {
                let entries = maplit::btreemap! {
                    ROOT_DIR  => (FileType::Directory, ".".to_owned()),
                    45        => (FileType::Directory, "..".to_owned()), // TODO
                    FASTA_DIR => (FileType::Directory, "fasta".to_owned()),
                    SEQ_DIR   => (FileType::Directory, "seqs".to_owned()),
                };
                for (o, (ino, entry)) in entries.iter().enumerate().skip(offset as usize) {
                    reply.add(*ino, o as i64 + 1, entry.0, entry.1.to_owned());
                }
                reply.ok();
            },
            FASTA_DIR => {
                let mut entries = maplit::btreemap! {
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
                let mut entries = maplit::btreemap! {
                    SEQ_DIR  => (FileType::Directory, ".".to_owned()),
                    ROOT_DIR => (FileType::Directory, "..".to_owned()),
                };
                for f in self.fragments.iter() {
                    entries.insert(f.seq_file.ino, (FileType::RegularFile, f.seq_file.name.clone()));
                }

                for (o, (ino, entry)) in entries.iter().enumerate().skip(offset as usize) {
                    debug!("{} {} {:?}", ino, o, entry);
                    reply.add(*ino, o as i64 + 1, entry.0, entry.1.to_owned());
                }
                reply.ok();
            },
            _ => {
                warn!("{} is not a directory", ino);
                reply.error(ENOENT);
            }
        }

    }

    fn unlink(&mut self, _req: &Request, parent: u64, name: &OsStr, reply: ReplyEmpty) {
        match parent {
            ROOT_DIR | SEQ_DIR | FASTA_DIR => {
                let name = name.to_str().unwrap();
                if self
                    .fragment_from_filename(name)
                    .and_then(|f| f.file_from_filename(name))
                    .is_some() {
                        self.fragments.retain(|f| f.fasta_file.name != name && f.seq_file.name != name);
                        reply.ok();
                    } else {
                        warn!("UNLINK: unknown file: `{:?}`", name);
                        reply.error(ENOENT);
                    }
            },
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
            // ADD_DIR = {
            // // TODO
            // // 1. Check if FASTA
            // // 2. Parse FASTA
            // // 3. Check for conflicts
            // // 4. Add new fragments to current fragments

            // let name = name.to_str().unwrap();
            // // If pathname already exists [...], this call fails with an EEXIST error.
            // if self.fragment_from_filename(name).is_some() {
            //     warn!("Cannot create `{:?}`, already exists", name);
            //     reply.error(EEXIST);
            //     return
            // }

            // let new_fragment = Fragment::with_empty_buffer(
            //     name, self.new_ino(), self.new_ino(),
            //     self.metadata.accessed().unwrap(), self.metadata.modified().unwrap(),
            // );

            // reply.entry(&TTL, &new_fragment.fasta_file.attrs, 0);
            // self.fragments.push(new_fragment);
            // },
            _ => {
                warn!("MKNOD: parent {} does not exist", parent);
                reply.error(ENOENT);
            }
        }
    }

    fn write(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, data: &[u8], _flags: u32, reply: ReplyWrite) {
        if self.fragment_from_ino(ino).is_none() {
            reply.error(ENOENT);
            return;
        }

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
    }

    fn setattr(&mut self, _req: &Request<'_>, ino: u64,
               mode: Option<u32>, uid: Option<u32>, gid: Option<u32>,
               size: Option<u64>,
               atime: Option<SystemTime>, mtime: Option<SystemTime>,
               _fh: Option<u64>,
               crtime: Option<SystemTime>, chgtime: Option<SystemTime>, bkuptime: Option<SystemTime>,
               flags: Option<u32>,
               reply: ReplyAttr) {
        info!("SETATTR");
        debug!("mode       {:?}", mode);
        debug!("gid        {:?}", gid);
        debug!("uid        {:?}", uid);
        debug!("size       {:?}", size);
        debug!("atime      {:?}", atime);
        debug!("mtime      {:?}", mtime);
        debug!("bkuptime   {:?}", bkuptime);
        debug!("chgtime    {:?}", chgtime);
        debug!("crtime     {:?}", crtime);
        debug!("flags      {:?}", flags);

        match ino {
            ROOT_DIR  => { reply.error(EACCES) }
            SEQ_DIR   => { reply.error(EACCES) }
            FASTA_DIR => { reply.error(EACCES) }
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
                                self.mut_fragment_from_ino(ino).map(|f| f.data = Backing::Buffer(Vec::new()));
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
                } else {
                    warn!("\t{:?} does not exist", ino);
                    reply.error(ENOENT);
                }
            }
        }
    }


    fn destroy(&mut self, _req: &Request) {}

    /// Rename a file.
    fn rename(&mut self, _req: &Request, parent: u64, name: &OsStr, newparent: u64, newname: &OsStr, reply: ReplyEmpty) {
        match parent {
            ROOT_DIR => {
                warn!("RENAME: access forbidden to ROOT_DIR");
                reply.error(EACCES);
            }
            SEQ_DIR | FASTA_DIR => {
                if newparent != parent {
                    error!("RENAME: Cannot move files oout of folder, please copy them");
                    reply.error(EACCES);
                } else {
                    if let Some(ref mut fragment) = self.mut_fragment_from_filename(name.to_str().unwrap()) {
                        fragment.id = newname.to_str().unwrap().to_string();
                        info!("RENAME: {:?} -> {:?}", name, newname);
                        self.concretize(true);
                        reply.ok()
                    } else {
                        warn!("\t{:?} does not exist", name);
                        reply.error(ENOENT);
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
        info!("FSYNC");
        self.concretize(false);
        reply.ok();
    }

    fn fsyncdir (&mut self, _req: &Request, _ino: u64, _fh: u64, _datasync: bool, reply: ReplyEmpty) {
        info!("FSYNCDIR");
        self.concretize(true);
        reply.ok();
    }

    fn flush(&mut self, _req: &Request, _ino: u64, _fh: u64, _lock_owner: u64, reply: ReplyEmpty) {
        info!("FLUSH");
        self.concretize(false);
        reply.ok();
    }

    /// Get file system statistics.
    fn statfs(&mut self, _req: &Request, _ino: u64, reply: ReplyStatfs) {
        info!("STATFS");
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

fn main() {
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
        .arg(Arg::with_name("DEST")
             .help("The mountpoint")
             .required(true)
             .index(2))

        .arg(Arg::with_name("v")
             .short("v")
             .multiple(true)
             .help("Sets the level of verbosity"))
        .arg(Arg::with_name("daemon")
             .short("D")
             .long("daemon")
             .help("Launch in the background; will exit when unmounted"))

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
        3 | _ => LevelFilter::Trace,
    };
    let log_config = ConfigBuilder::new()
        .set_time_format_str("")
        .build();

    if args.is_present("daemon") {
        let mut log_file = std::env::temp_dir();
        log_file.push("fusta.log");
        WriteLogger::new(log_level, log_config, std::fs::File::create(log_file).unwrap());

        let mut pid_file = std::env::temp_dir();
        pid_file.push("fusta.pid");

        Daemonize::new()
            .pid_file(pid_file)
            .working_directory(std::env::current_dir().expect("Unable to read current directory"))
            .start().unwrap();
    } else {
        TermLogger::init(log_level, log_config, TerminalMode::Mixed).expect("Unable to initialize logger");
    }


    let fasta_file = value_t!(args, "FASTA", String).unwrap();
    let mountpoint = value_t!(args, "DEST", String).unwrap();
    let mut fuse_options: Vec<&OsStr> = vec![
        &OsStr::new("-o"), &OsStr::new("auto_unmount"),
        &OsStr::new("-o"), &OsStr::new("atomic_o_trunc"),
        &OsStr::new("-o"), &OsStr::new("default_permissions"),
        &OsStr::new("-o"), &OsStr::new("subtype"), &OsStr::new(&fasta_file),
    ];
    let settings = FustaSettings {
        mmap: args.is_present("mmap"),
    };

    info!("Using MMAP:      {}", settings.mmap);

    let fs = FustaFS::new(settings, &fasta_file);
    if !std::path::Path::new(&mountpoint).is_dir() {
        error!("{} is not a directory", mountpoint);
        std::process::exit(1);
    }
    if args.is_present("nonempty") {
        fuse_options.push(&OsStr::new("-o"));
        fuse_options.push(&OsStr::new("nonempty"));
    }
    if !args.is_present("nonempty") && std::fs::read_dir(&mountpoint).unwrap().take(1).count() != 0 {
        error!("{} is not empty. Use the -E flag if you want to mount in a non-empty directory.", mountpoint);
        std::process::exit(1);
    }

    ctrlc::set_handler(move || {
        info!("Ctrl-C received, exiting.");
        std::process::exit(0);
    }).expect("Error setting Ctrl-C handler");

    fuse::mount(fs, &mountpoint, &fuse_options).unwrap();
}
