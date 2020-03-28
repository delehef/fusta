use std::env;
use fuse::{FileType, FileAttr, Filesystem, Request, ReplyData, ReplyEntry, ReplyAttr, ReplyDirectory};

use std::ffi::OsStr;
use std::time::Duration;
use std::fs;
use libc::ENOENT;

use memmap;
use std::io::SeekFrom;
use std::io::prelude::*;

use clap::*;

use fusta::fasta::*;

const TTL: Duration = Duration::from_secs(100);


struct FustaSettings {
    mmap: bool,
    labels: bool,
}

struct FustaFS {
    fasta: Vec<Fragment>,
    metadata: fs::Metadata,
    root_dir_attrs: FileAttr,
    entries: Vec<(u64, fuse::FileType, String)>,
    filename: String,
    mmap: Option<memmap::Mmap>,
    settings: FustaSettings,
}


impl FustaFS {
    fn new(settings: FustaSettings, filename: &str) -> FustaFS {
        eprintln!("Reading {}...", filename);
        let fasta = from_file(filename).unwrap();
        eprintln!("Done.");
        let mut keys = fasta.iter().map(|f| &f.id).collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        if keys.len() != fasta.len() {panic!("Duplicated keys")}

        let mut entries = vec![
            (1, FileType::Directory, ".".to_owned()),
            (1, FileType::Directory, "..".to_owned()),
        ];
        for (i, fragment) in fasta.iter().enumerate() {
            entries.push(((i + 2) as u64, FileType::RegularFile, fragment.id.to_owned()))
        }

        let f = std::fs::File::open(filename).unwrap();

        FustaFS {
            fasta: fasta,
            filename: filename.to_owned(),
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
            mmap: if settings.mmap { Some(unsafe { memmap::MmapOptions::new().map(&f).unwrap() })} else { None },
            settings: settings,
        }
    }

    fn fragment_to_fileattrs(&self, ino: u64, fragment: &Fragment) -> FileAttr {
        FileAttr {
            ino: ino,
            size: (fragment.pos.1 - fragment.pos.0) as u64,
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
        eprintln!("LOOKUP {}/{:?}", parent, name);

        if parent == 1  {
            if let Some(i) = (0 .. self.fasta.len()).find(|&i| self.fasta[i].id == name.to_str().unwrap()) {
                let fragment = &self.fasta[i];
                let ino = i as u64 + 2;
                reply.entry(&TTL, &self.fragment_to_fileattrs(ino, &fragment), 0);
            } else {
                eprintln!("\t{:?} does not exist", name);
                reply.error(ENOENT);
            }
        } else {
            eprintln!("\tParent {} does not exist", parent);
            reply.error(ENOENT);
        }
    }

    fn getattr(&mut self, _req: &Request, ino: u64, reply: ReplyAttr) {
        eprintln!("GETATTR {}", ino);

        match ino {
            1 => reply.attr(&TTL, &self.root_dir_attrs),
            ino if ino >= 2 && ino < self.fasta.len() as u64 + 2 => {
                let fasta_id = (ino - 2) as usize;
                let fragment = &self.fasta[fasta_id];
                reply.attr(&TTL, &self.fragment_to_fileattrs(ino, &fragment))
            },
            _ => {
                eprintln!("\tino {} does not exist ({:?})", ino, self.fasta.len());
                reply.error(ENOENT)
            },
        }
    }

    fn read(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, size: u32, reply: ReplyData) {
        eprintln!("READ {} {} | {}", ino, offset, size);

        if ino >= 2 && ino <= self.fasta.len() as u64 + 2 {
            let id = (ino - 2) as usize;

            let mut buffer;
            let out = if let Some(mmap) = &self.mmap {
                eprintln!("Using MMAP");

                let end = if offset + size as i64 >= self.fasta[id].len as i64 {
                    self.fasta[id].len as i64
                } else {
                    offset + size as i64
                };
                &mmap[self.fasta[id].pos.0 + offset as usize .. self.fasta[id].pos.0 + end as usize]
            } else { // fseek
                buffer = vec![0u8; size as usize];
                let mut f = std::fs::File::open(&self.filename).unwrap();
                f.seek(SeekFrom::Start(self.fasta[id].pos.0 as u64 + offset as u64)).unwrap();
                let _ = f.read(&mut buffer).unwrap();
                &buffer
            };

            reply.data(&out)

        } else {
            eprintln!("\t{} is not a file", ino);
            reply.error(ENOENT);
        }
    }

    fn readdir(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, mut reply: ReplyDirectory) {
        eprintln!("READDIR {}", ino);

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
        .arg(Arg::with_name("mmap")
             .short("M")
             .long("mmap")
             .help("Use mmap instead of seek to browse FASTA fragments. Faster, but memory hungry"))

        .arg(Arg::with_name("nonempty")
             .short("E")
             .long("non-empty")
             .help("Perform the mount even if the destination folder is not empty"))

        .get_matches();

    let fasta_file = value_t!(args, "FASTA", String).unwrap();
    let mountpoint = value_t!(args, "DEST", String).unwrap();
    let mut fuse_options: Vec<&OsStr> = Vec::new();
    let settings = FustaSettings {
        mmap: args.is_present("mmap"),
        labels: false,
    };


    let fs = FustaFS::new(settings, &fasta_file);
    if !std::path::Path::new(&mountpoint).is_dir() {
        panic!("{} is not a directory", mountpoint)
    }

    if args.is_present("nonempty") {
        fuse_options.push(&OsStr::new("-o"));
        fuse_options.push(&OsStr::new("nonempty"));
    }
    if !args.is_present("nonempty") && std::fs::read_dir(&mountpoint).unwrap().take(1).count() != 0 {
        panic!("{} is not empty. Use the -E flag if you want to mount in a non-empty directory.", mountpoint);
    }

    eprintln!("Options => {:?}", fuse_options);

    fuse::mount(fs, &mountpoint, &fuse_options).unwrap();
}
