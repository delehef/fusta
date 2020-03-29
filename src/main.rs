use std::env;
use fuse::{FileType, FileAttr, Filesystem, Request, ReplyData, ReplyEntry, ReplyAttr, ReplyDirectory};

use std::ffi::OsStr;
use std::time::Duration;
use std::fs;
use libc::ENOENT;
use daemonize::*;

use memmap;
use std::io::SeekFrom;
use std::io::prelude::*;

use clap::*;

use fusta::fasta::*;

use log::*;
use simplelog::*;

const TTL: Duration = Duration::from_secs(100);
const INO_DIR: u64 = 1;

#[derive(Debug)]
enum Backing {
    File(String, usize, usize),
    Buffer(Vec<u8>),
    MMap(memmap::Mmap)
}

struct BufferedFragment {
    id: String,
    name: Option<String>,
    data: Backing,
}
impl BufferedFragment {
    fn len(&self) -> usize {
        match &self.data {
            Backing::File(_, start, end) => end - start,
            Backing::Buffer(ref b)       => b.len(),
            Backing::MMap(ref mmap)          => mmap.len()
        }
    }

    fn data(&self) -> Box<[u8]> {
        match &self.data {
            Backing::File(filename, start, _) => {
                let mut buffer = vec![0u8; self.len() as usize];
                let mut f = std::fs::File::open(&filename).unwrap();
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
                let mut f = std::fs::File::open(&filename).unwrap();
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
}

struct FustaSettings {
    mmap: bool,
    keep_labels: bool,
}

struct FustaFS {
    fragments: Vec<BufferedFragment>,
    metadata: fs::Metadata,
    root_dir_attrs: FileAttr,
    entries: Vec<(u64, fuse::FileType, String)>,
    filename: String,
    settings: FustaSettings,
}


impl FustaFS {
    fn new(settings: FustaSettings, filename: &str) -> FustaFS {
        info!("Reading {}...", filename);
        let fasta_file = std::fs::File::open(filename).unwrap();
        let fragments = FastaReader::new(fasta_file, settings.keep_labels).collect::<Vec<_>>();
        info!("Done.");
        let mut keys = fragments.iter().map(|f| &f.id).collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        if keys.len() != fragments.len() {panic!("Duplicated keys")}

        let mut entries = vec![
            (INO_DIR, FileType::Directory, ".".to_owned()),
            (INO_DIR, FileType::Directory, "..".to_owned()),
        ];
        for (i, fragment) in fragments.iter().enumerate() {
            entries.push(((i + 2) as u64, FileType::RegularFile, fragment.id.to_owned()))
        }

        let file = std::fs::File::open(filename).unwrap();

        FustaFS {
            fragments: fragments.iter().map(|f| BufferedFragment {
                id: f.id.clone(),
                name: f.name.clone(),
                data: if settings.mmap {
                    Backing::MMap(unsafe { memmap::MmapOptions::new().offset(f.pos.0 as u64).len(f.len).map(&file).unwrap() })
                } else {
                    Backing::File(filename.to_owned(), f.pos.0, f.pos.1)
                }

            }).collect::<Vec<_>>(),
            filename: filename.to_owned(),
            root_dir_attrs: FileAttr {
                ino: INO_DIR,
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
            settings: settings,
        }
    }

    fn fragment_to_fileattrs(&self, ino: u64, fragment: &BufferedFragment) -> FileAttr {
        FileAttr {
            ino: ino,
            size: fragment.len() as u64,
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
        debug!("LOOKUP {}/{:?}", parent, name);

        if parent == INO_DIR  {
            if let Some(i) = (0 .. self.fragments.len()).find(|&i| self.fragments[i].id == name.to_str().unwrap()) {
                let fragment = &self.fragments[i];
                let ino = i as u64 + 2;
                reply.entry(&TTL, &self.fragment_to_fileattrs(ino, &fragment), 0);
            } else {
                debug!("\t{:?} does not exist", name);
                reply.error(ENOENT);
            }
        } else {
            debug!("\tParent {} does not exist", parent);
            reply.error(ENOENT);
        }
    }

    fn getattr(&mut self, _req: &Request, ino: u64, reply: ReplyAttr) {
        debug!("GETATTR {}", ino);

        match ino {
            INO_DIR => reply.attr(&TTL, &self.root_dir_attrs),
            ino if ino >= 2 && ino < self.fragments.len() as u64 + 2 => {
                let fasta_id = (ino - 2) as usize;
                let fragment = &self.fragments[fasta_id];
                reply.attr(&TTL, &self.fragment_to_fileattrs(ino, &fragment))
            },
            _ => {
                debug!("\tino {} does not exist ({:?})", ino, self.fragments.len());
                reply.error(ENOENT)
            },
        }
    }

    fn read(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, size: u32, reply: ReplyData) {
        debug!("READ {} {} | {}", ino, offset, size);

        if ino >= 2 && ino <= self.fragments.len() as u64 + 2 {
            let id = (ino - 2) as usize;
            reply.data(&self.fragments[id].chunk(offset, size))
        } else {
            debug!("\t{} is not a file", ino);
            reply.error(ENOENT);
        }
    }

    fn readdir(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, mut reply: ReplyDirectory) {
        debug!("READDIR {}", ino);

        if ino != INO_DIR {
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
        .arg(Arg::with_name("labels")
             .short("L")
             .long("keep-labels")
             .help("Keep the FASTA labels in the virtual files"))
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
        let mut log_file = env::temp_dir();
        log_file.push("fusta.log");
        WriteLogger::new(log_level, log_config, std::fs::File::create(log_file).unwrap());

        let mut pid_file = env::temp_dir();
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
    let mut fuse_options: Vec<&OsStr> = vec![&OsStr::new("-o"), &OsStr::new("auto_unmount")];
    let settings = FustaSettings {
        mmap: args.is_present("mmap"),
        keep_labels: args.is_present("labels"),
    };

    info!("Using MMAP:      {}", settings.mmap);
    info!("Keeping labels:  {}", settings.keep_labels);

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
