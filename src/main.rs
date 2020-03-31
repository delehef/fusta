use std::env;
use fuse::*;

use std::ffi::OsStr;
use std::time::{SystemTime, Duration};
use std::fs;
use libc::*;
use daemonize::*;
use std::collections::BTreeMap;

use memmap;
use std::io::SeekFrom;
use std::io::prelude::*;

use clap::*;

use fusta::fasta::*;

use log::*;
use simplelog::*;

const TTL: Duration = Duration::from_secs(1);
const INO_DIR: u64 = 1;
const FASTA_EXT: &str = "fa";
const EXTENSIONS: &[&str] = &["FA", "FASTA"];

#[derive(Debug)]
enum Backing {
    File(String, usize, usize),
    Buffer(Vec<u8>),
    MMap(memmap::Mmap)
}

#[derive(Debug)]
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
}

struct FustaSettings {
    mmap: bool,
    keep_labels: bool,
}

struct FustaFS {
    fragments: BTreeMap<u64, BufferedFragment>,
    metadata: fs::Metadata,
    root_dir_attrs: FileAttr,
    files_shared_attrs: FileAttr,
    filename: String,
    settings: FustaSettings,
}


impl FustaFS {
    fn new(settings: FustaSettings, filename: &str) -> FustaFS {
        let metadata = fs::metadata(filename).unwrap();
        let mut r = FustaFS {
            fragments: BTreeMap::new(),
            filename: String::new(),
            root_dir_attrs: FileAttr {
                ino: INO_DIR,
                size: 0,
                blocks: 0,
                atime: std::time::SystemTime::now(),
                mtime: std::time::SystemTime::now(),
                ctime: std::time::SystemTime::now(),
                crtime: std::time::SystemTime::now(),
                kind: FileType::Directory,
                perm: 0o775,
                nlink: 1,
                uid: unsafe { libc::geteuid() },
                gid: unsafe { libc::getgid() },
                rdev: 0,
                flags: 0,
            },
            files_shared_attrs: FileAttr {
                ino: 0,
                size: 0,
                blocks: 1,
                atime: metadata.accessed().unwrap(),
                mtime: metadata.modified().unwrap(),
                ctime: metadata.modified().unwrap(),
                crtime: metadata.modified().unwrap(),
                kind: FileType::RegularFile,
                perm: 0o664,
                nlink: 1,
                uid: unsafe { libc::geteuid() },
                gid: 20,
                rdev: 0,
                flags: 0,
            },
            metadata: metadata,
            settings: settings,
        };

        r.read_fasta(filename);
        r
    }

    fn new_ino(&self) -> u64 {
        let r =  if self.fragments.keys().next().is_none() {
            3
        } else {
            self.fragments.keys().max_by(|k1, k2| k1.cmp(k2)).unwrap() + 1
        };
        debug!("New ino: {}", r);
        r
    }

    fn read_fasta(&mut self, filename: &str) {
        info!("Reading {}...", filename);
        let fasta_file = fs::File::open(filename).unwrap();
        let fragments = FastaReader::new(fasta_file, self.settings.keep_labels).collect::<Vec<_>>();
        info!("Done.");
        let mut keys = fragments.iter().map(|f| &f.id).collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        if keys.len() != fragments.len() {panic!("Duplicated keys")}

        let file = fs::File::open(filename).unwrap();
        self.filename = filename.to_owned();

        self.fragments.clear();
        for f in fragments.iter() {
            self.fragments.insert(
                self.new_ino(),
                BufferedFragment {
                    id: f.id.clone(),
                    name: f.name.clone(),
                    data: if self.settings.mmap {
                        Backing::MMap(unsafe { memmap::MmapOptions::new().offset(f.pos.0 as u64).len(f.len).map(&file).unwrap() })
                    } else {
                        Backing::File(filename.to_owned(), f.pos.0, f.pos.1)
                    }
                }
            );
        }
    }

    fn concretize(&mut self, force: bool) {
        // if force is set to true, we still force the rewrite so that "phantom" (deleted) fragments
        // are effectively removed
        if !force && !self.fragments.values().any(|f| matches!(f.data, Backing::Buffer(_))) {
            debug!("Nothing to concretize; leaving");
            return;
        }

        let tmp_filename = format!("{}#fusta#", &self.filename);
        { // Scope to ensure the content is correctly written
            trace!("Writing fragments");
            let mut tmp_file = fs::File::create(&tmp_filename).unwrap();
            for fragment in self.fragments.values() {
                trace!("Writing {}", fragment.id);
                if !self.settings.keep_labels {
                    tmp_file.write_all(format!(
                        ">{} {}",
                        fragment.id, fragment.name.as_ref().unwrap_or(&String::new())
                    ).as_bytes()).unwrap();
                }

                tmp_file.write_all(&fragment.data()).unwrap();
            }
        }
        trace!("Renaming {} to {}", tmp_filename, &self.filename);
        fs::rename(tmp_filename, &self.filename).unwrap();
        trace!("Rebuilding index");
        let filename = self.filename.clone();
        self.read_fasta(&filename);
    }

    fn fragment_to_fileattrs(&self, ino: u64, fragment: &BufferedFragment) -> FileAttr {
        FileAttr {
            ino: ino,
            size: fragment.len() as u64,
            .. self.files_shared_attrs
        }
    }

    fn fragment_from_ino(&mut self, ino: u64) -> Option<&mut BufferedFragment> {
        self.fragments.get_mut(&ino)
    }

    fn fragment_from_filename(&self, name: &OsStr) -> Option<(&u64, &BufferedFragment)> {
        self.fragments.iter().find(|(_, f)| f.to_filename() == name.to_str().unwrap())
    }

    fn mut_fragment_from_filename(&mut self, name: &OsStr) -> Option<(&u64, &mut BufferedFragment)> {
        self.fragments.iter_mut().find(|(_, f)| f.to_filename() == name.to_str().unwrap())
    }
}


impl Filesystem for FustaFS {
    fn lookup(&mut self, _req: &Request, parent: u64, name: &OsStr, reply: ReplyEntry) {
        debug!("LOOKUP {}/{:?}", parent, name);
        if parent != INO_DIR  {
            warn!("\tParent {} does not exist", parent);
            reply.error(ENOENT);
            return;
        }

        if let Some((&i, fragment)) = self.fragment_from_filename(name) {
            reply.entry(&TTL, &self.fragment_to_fileattrs(i, fragment), 0);
        } else {
            warn!("\t{:?} does not exist", name);
            reply.error(ENOENT);
        }
    }

    fn getattr(&mut self, _req: &Request, ino: u64, reply: ReplyAttr) {
        debug!("GETATTR {}", ino);

        match ino {
            INO_DIR => reply.attr(&TTL, &self.root_dir_attrs),
            ino if self.fragments.contains_key(&ino) => {
                let fragment = &self.fragments[&ino];
                reply.attr(&TTL, &self.fragment_to_fileattrs(ino, &fragment))
            },
            _ => {
                warn!("\tino {} does not exist ({:?})", ino, self.fragments.len());
                reply.error(ENOENT)
            },
        }
    }

    fn read(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, size: u32, reply: ReplyData) {
        debug!("READ {} {} | {}", ino, offset, size);

        if self.fragments.contains_key(&ino) {
            reply.data(&self.fragments[&ino].chunk(offset, size))
        } else {
            warn!("\t{} is not a file", ino);
            reply.error(ENOENT);
        }
    }

    fn readdir(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, mut reply: ReplyDirectory) {
        debug!("READDIR {}", ino);

        if ino != INO_DIR {
            warn!("{} is not a directory", ino);
            reply.error(ENOENT);
            return;
        }

        let mut entries = maplit::btreemap! {
            INO_DIR     => (FileType::Directory, ".".to_owned()),
            INO_DIR + 1 => (FileType::Directory, "..".to_owned()),
        };
        for (&i, fragment) in self.fragments.iter() {
            entries.insert(i, (FileType::RegularFile, fragment.to_filename()));
        }

        for (o, (ino, entry)) in entries.iter().enumerate().skip(offset as usize) {
            reply.add(*ino, o as i64 + 1, entry.0, entry.1.to_owned());
        }
        reply.ok();
    }

    fn unlink(&mut self, _req: &Request, parent: u64, name: &OsStr, reply: ReplyEmpty) {
        debug!("UNLINK {}/{:?}", parent, name);
        if parent != INO_DIR {
            warn!("{} is not a directory", parent);
            reply.error(ENOENT);
            return;
        }

        if let Some((&ino, _)) = self.fragment_from_filename(name) {
            match self.fragments.remove(&ino) {
                Some(_) => {
                    self.concretize(true);
                    reply.ok();
                },
                None    => {
                    warn!("Unknown ino: `{:?}`", name);
                    reply.error(ENOENT);
                }
            }
        } else {
            warn!("Unknown file: `{:?}`", name);
            reply.error(ENOENT);
        }
    }

    fn mknod(&mut self, _req: &Request, parent: u64, name: &OsStr, _mode: u32, _rdev: u32, reply: ReplyEntry) {
        if parent != INO_DIR {
            reply.error(ENOENT);
            return;
        }

        // If pathname already exists [...], this call fails with an EEXIST error.
        if self.fragment_from_filename(name).is_some() {
            warn!("Cannot create `{:?}`, already exists", name);
            reply.error(EEXIST);
            return
        }

        let new_ino = self.new_ino();
        let new_fragment = BufferedFragment {
            id: name.to_str().unwrap().to_string(),
            name: None,
            data: Backing::Buffer(format!(">{:?}\n", name).as_bytes().to_vec())
        };
            reply.entry(&TTL, &self.fragment_to_fileattrs(new_ino, &new_fragment), 0);
        self.fragments.insert(new_ino, new_fragment);
    }

    fn create(&mut self, _req: &Request, _parent: u64, name: &OsStr, _mode: u32, flags: u32, reply: ReplyCreate) {
        warn!("CREATE");
        let name = name.to_str().unwrap().to_string();

        let new_ino = self.new_ino();
        let new_fragment = BufferedFragment {
            id: name,
            name: None,
            data: Backing::Buffer(Vec::new())
        };

        reply.created(&TTL, &self.fragment_to_fileattrs(new_ino, &new_fragment), 0, 0, flags);
        self.fragments.insert(new_ino, new_fragment);
    }

    fn write(&mut self, _req: &Request, ino: u64, _fh: u64, offset: i64, data: &[u8], _flags: u32, reply: ReplyWrite) {
        debug!("WRITE {}: {}/{:?}", ino, offset, data);
        if let Some(mut fragment) = self.fragment_from_ino(ino) {
            // As soon as there's a write, we have to switch to a buffer-backed storage
            if !matches!(fragment.data, Backing::Buffer(_)) {
                fragment.data = Backing::Buffer(fragment.data().to_vec());
            }

            // Ensure that the backing buffer is big enough
            let max_size = offset as usize + data.len();
            if max_size > fragment.len() { fragment.extend(max_size) }

            // Then finally write the data
            if let Backing::Buffer(b) = &mut fragment.data {
                let start = offset as usize;
                let end = start + data.len();
                b.splice(start .. end, data.iter().cloned());
                reply.written(data.len() as u32);
            } else {
                panic!("Something went very wrong...")
            }
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
        error!("SETATTR");
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

        if let Some(ref mut fragment) = self.fragments.get_mut(&ino) {
            if let Some(uid) = uid         { self.files_shared_attrs.uid = uid }
            if let Some(gid) = gid         { self.files_shared_attrs.gid = gid }
            // TODO
            // - store per file
            // - see how setuid is suid/sgid works
            if let Some(mode) = mode       { self.files_shared_attrs.perm = mode as u16 }
            if let Some(atime) = atime     { self.files_shared_attrs.atime = atime }
            if let Some(mtime) = mtime     { self.files_shared_attrs.mtime = mtime }
            if let Some(chgtime) = chgtime { self.files_shared_attrs.mtime = chgtime }
            if let Some(crtime) = crtime   { self.files_shared_attrs.crtime = crtime }
            // macOS only
            if let Some(flags) = flags     { self.files_shared_attrs.flags = flags }
            if let Some(size) = size {
                let size = size as usize;
                if size == 0 { // Clear the file, called by the truncate syscall
                    fragment.data = Backing::Buffer(Vec::new());
                } else if size  != fragment.len() { // Redim the file
                    if !matches!(fragment.data, Backing::Buffer(_)) {
                        fragment.data = Backing::Buffer(fragment.data().to_vec());
                    }
                    fragment.extend(size)
                }
                self.files_shared_attrs.size = size as u64;
            }

            reply.attr(&TTL, &self.fragment_to_fileattrs(ino, &self.fragments[&ino]));
        } else {
            warn!("\t{:?} does not exist", ino);
            reply.error(ENOENT);
        }
    }

    fn destroy(&mut self, _req: &Request) {}

    /// Rename a file.
    fn rename(&mut self, _req: &Request, _parent: u64, _name: &OsStr, _newparent: u64, _newname: &OsStr, reply: ReplyEmpty) {
        error!("RENAME");
        reply.error(ENOSYS);
    }

    fn fsync(&mut self, _req: &Request, _ino: u64, _fh: u64, _datasync: bool, reply: ReplyEmpty) {
        info!("FSYNC");
        self.concretize(false);
        reply.ok();
    }

    fn fsyncdir (&mut self, _req: &Request, _ino: u64, _fh: u64, _datasync: bool, reply: ReplyEmpty) {
        info!("FSYNCDIR");
        self.concretize(true);
        reply.error(ENOSYS);
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
    let mut fuse_options: Vec<&OsStr> = vec![
        &OsStr::new("-o"), &OsStr::new("auto_unmount"),
        &OsStr::new("-o"), &OsStr::new("atomic_o_trunc"),
        &OsStr::new("-o"), &OsStr::new("default_permissions"),
        &OsStr::new("-o"), &OsStr::new("subtype"), &OsStr::new(&fasta_file),
    ];
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
