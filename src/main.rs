use std::env;
use fuse::{FileType, FileAttr, Filesystem, Request, ReplyData, ReplyEntry, ReplyAttr, ReplyDirectory};

use std::ffi::OsStr;
use std::time::Duration;
use std::fs;
use libc::ENOENT;

use fusta::fasta::*;

const TTL: Duration = Duration::from_secs(100);


struct FustaFS {
    fasta: Vec<Fragment>,
    metadata: fs::Metadata,
    root_dir_attrs: FileAttr,
    entries: Vec<(u64, fuse::FileType, String)>,
}


impl FustaFS {
    fn new(filename: &str) -> FustaFS {
        let fasta = from_file(filename).unwrap();
        let mut entries = vec![
            (1, FileType::Directory, ".".to_owned()),
            (1, FileType::Directory, "..".to_owned()),
        ];
        for (i, fragment) in fasta.iter().enumerate() {
            entries.push(((i + 2) as u64, FileType::RegularFile, fragment.id.to_owned()))
        }

        FustaFS {
            fasta: fasta,
            root_dir_attrs: FileAttr {
                ino: 1,
                size: 0,
                blocks: 0,
                atime: std::time::SystemTime::now(),
                mtime: std::time::SystemTime::now(),
                ctime: std::time::SystemTime::now(),
                crtime: std::time::SystemTime::now(),
                kind: FileType::Directory,
                perm: 0o555,
                nlink: 2,
                uid: unsafe { libc::geteuid() },
                gid: unsafe { libc::getgid() },
                rdev: 0,
                flags: 0,
            },
            metadata: fs::metadata(filename).unwrap(),
            entries: entries,
        }
    }

    fn fragment_to_fileattrs(&self, ino: u64, fragment: &Fragment) -> FileAttr {
        FileAttr {
            ino: ino,
            size: fragment.sequence.len() as u64,
            blocks: 1,
            atime: self.metadata.accessed().unwrap(),
            mtime: self.metadata.modified().unwrap(),
            ctime: self.metadata.modified().unwrap(),
            crtime: self.metadata.modified().unwrap(),
            kind: FileType::RegularFile,
            perm: 0o444,
            nlink: 1,
            uid: unsafe { libc::geteuid() },
            gid: 20,
            rdev: 0,
            flags: 0,
        }
    }
}


impl Filesystem for FustaFS {
    fn lookup(&mut self, _req: &Request, parent: u64, name: &OsStr, reply: ReplyEntry) {
        if parent == 1  {
            if let Some(i) = (0 .. self.fasta.len()).find(|&i| self.fasta[i].id == name.to_str().unwrap()) {
                let fragment = &self.fasta[i];
                let ino = i as u64 + 2;
                reply.entry(&TTL, &self.fragment_to_fileattrs(ino, &fragment), 0);
            } else {
                reply.error(ENOENT);
            }
        } else {
            reply.error(ENOENT);
        }
    }

    fn getattr(&mut self, _req: &Request, ino: u64, reply: ReplyAttr) {
        match ino {
            1 => reply.attr(&TTL, &self.root_dir_attrs),
            ino if ino >= 2 && ino < self.fasta.len() as u64 - 2 => {
                let fasta_id = (ino - 2) as usize;
                let fragment = &self.fasta[fasta_id];
                reply.attr(&TTL, &self.fragment_to_fileattrs(ino, &fragment))
            },
            _ => reply.error(ENOENT),
        }
    }

    fn read(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, size: u32, reply: ReplyData) {
        if ino >= 2 && ino <= self.fasta.len() as u64 + 2 {
            let id = (ino - 2) as usize;
            eprintln!("{} -> {}/{}B\n{} & {}", ino, id, self.fasta[id].sequence.len(), offset, size);
            let end = std::cmp::min(self.fasta[id].sequence.len() as u32, offset as u32 + size) as usize;
            reply.data(&self.fasta[id].sequence[offset as usize .. end]);
        } else {
            reply.error(ENOENT);
        }
    }

    fn readdir(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, mut reply: ReplyDirectory) {
        if ino != 1 {
            reply.error(ENOENT);
            return;
        }

        for (i, entry) in self.entries.iter().enumerate().skip(offset as usize) {
            reply.add(entry.0, (i+1) as i64, entry.1, entry.2.to_owned());
        }
        reply.ok();
    }
}

fn main() {
    env_logger::init();
    let mountpoint = env::args_os().nth(1).unwrap();
    let filename = &env::args().nth(2).unwrap();
    eprintln!("Reading {}...", filename);
    let fs = FustaFS::new(filename);
    eprintln!("Done.");

    fuse::mount(fs, &mountpoint, &[]).unwrap();
}
