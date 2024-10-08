#![allow(clippy::redundant_field_names)]
use crate::notify;
use anyhow::{Context, Result};
use fuser::*;
use libc::*;
use log::*;
use maplit::*;
use multi_map::MultiMap;
use regex::Regex;
use smartstring::SmartString;
use std::cell::RefCell;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::convert::TryInto;
use std::ffi::OsStr;
use std::fs;
use std::time::{Duration, SystemTime};
type SString = SmartString<smartstring::LazyCompact>;

use std::io::prelude::*;
use std::io::SeekFrom;

use fusta::fasta::*;

#[cfg(unix)]
const FORBIDDEN_CHARS: [char; 2] = ['\\', '\0'];
#[cfg(windows)]
const FORBIDDEN_CHARS: [char; 8] = ['\\', '/', ':', '*', '?', '|', '<', '>'];
const TTL: Duration = Duration::from_secs(1);

const FASTA_EXT: &str = ".fa";
const SEQ_EXT: &str = ".seq";
const DOTLESS_SEQ_EXT: &str = ".seq";

// Virtual directories
const ROOT_DIR: u64 = 1;
const FASTA_DIR: u64 = 2;
const SEQ_DIR: u64 = 3;
const APPEND_DIR: u64 = 4;
const SUBFRAGMENTS_DIR: u64 = 5;

// First free ino
const FIRST_INO: u64 = 20;

// Pure virtual files
const INFO_FILE: u64 = 10;
const INFO_FILE_NAME: &str = "infos.txt";
const LABELS_FILE: u64 = 11;
const LABELS_FILE_NAME: &str = "labels.txt";
const INFO_CSV_FILE: u64 = 12;
const INFO_CSV_FILE_NAME: &str = "infos.csv";

fn is_fasta_char(c: u8) -> bool {
    c.is_ascii_alphanumeric() || [b'\n', b'-', b'_', b'.', b'+', b'='].contains(&c)
}

#[derive(Debug, PartialEq, Clone)]
enum FileClass {
    Fasta(std::cell::RefCell<Vec<u8>>),
    Seq,
    Text,
}

trait VirtualFile {
    fn name(&self) -> &str;
    fn attrs(&self) -> &FileAttr;
    fn mut_attrs(&mut self) -> &mut FileAttr;
    fn ino(&self) -> u64;
    fn class(&self) -> &FileClass;
    fn data(&self) -> &[u8];
    fn set_data(&mut self, data: &[u8]);
}

struct BufferFile {
    name: SString,
    ino: u64,
    attrs: FileAttr,
    class: FileClass,
    _data: Vec<u8>,
}
impl VirtualFile for BufferFile {
    fn name(&self) -> &str {
        &self.name
    }
    fn attrs(&self) -> &FileAttr {
        &self.attrs
    }
    fn mut_attrs(&mut self) -> &mut FileAttr {
        &mut self.attrs
    }
    fn ino(&self) -> u64 {
        self.ino
    }
    fn class(&self) -> &FileClass {
        &self.class
    }
    fn data(&self) -> &[u8] {
        &self._data
    }
    fn set_data(&mut self, data: &[u8]) {
        self._data.clear();
        self._data.extend_from_slice(data)
    }
}

#[derive(Debug)]
struct FragmentFile {
    name: SString,
    ino: u64,
    attrs: FileAttr,
    class: FileClass,
}
impl VirtualFile for FragmentFile {
    fn name(&self) -> &str {
        &self.name
    }
    fn attrs(&self) -> &FileAttr {
        &self.attrs
    }
    fn mut_attrs(&mut self) -> &mut FileAttr {
        &mut self.attrs
    }
    fn ino(&self) -> u64 {
        self.ino
    }
    fn class(&self) -> &FileClass {
        &self.class
    }
    fn data(&self) -> &[u8] {
        unimplemented!()
    }
    fn set_data(&mut self, _: &[u8]) {
        unimplemented!()
    }
}

#[derive(Debug)]
enum Backing {
    File(SString, usize, usize), // A start, end pair in a file
    Buffer(Vec<u8>),             // A chunk of memory
    PureBuffer(Vec<u8>), // The same, but guaranteed pure (i.e. no newlines) - can be accessed directly
    MMap(memmap2::Mmap), // A memmapped chunk of memory
}
impl Backing {
    fn len(&self) -> usize {
        match self {
            Backing::File(_, start, end) => end - start,
            Backing::Buffer(ref b) => b.len(),
            Backing::PureBuffer(ref b) => b.len(),
            Backing::MMap(ref mmap) => mmap.len(),
        }
    }
}

#[derive(Debug)]
struct Fragment {
    id: SString,
    name: Option<String>,
    data: Backing,
    fasta_file: FragmentFile,
    seq_file: FragmentFile,
}
impl Fragment {
    fn make_label(id: &str, name: &Option<String>) -> String {
        format!(
            ">{}{}\n",
            id,
            name.as_ref().map(|n| format!(" {}", n)).unwrap_or_default()
        )
    }

    fn make_virtual_file(
        ino: u64,
        name: &str,
        permissions: u16,
        size: usize,
        class: FileClass,
        accessed: SystemTime,
        modified: SystemTime,
    ) -> FragmentFile {
        FragmentFile {
            name: name.into(),
            ino: ino,
            class: class,
            attrs: FileAttr {
                ino: ino,
                size: size as u64,
                blocks: 1,
                atime: accessed,
                mtime: modified,
                ctime: modified,
                crtime: modified,
                kind: FileType::RegularFile,
                perm: permissions,
                nlink: 0,
                uid: unsafe { libc::geteuid() },
                gid: unsafe { libc::getgid() },
                rdev: 0,
                flags: 0,
                blksize: 512,
            },
        }
    }

    fn new(
        id: &str,
        name: &Option<String>,
        data: Backing,
        fasta_ino: u64,
        seq_ino: u64,
        accessed: SystemTime,
        modified: SystemTime,
    ) -> Fragment {
        let label = Fragment::make_label(id, name);
        let data_size = data.len();
        Fragment {
            id: id.into(),
            name: name.clone(),
            data: data,
            fasta_file: Fragment::make_virtual_file(
                fasta_ino,
                &format!("{}{}", id, FASTA_EXT),
                0o444,
                label.as_bytes().len() + data_size,
                FileClass::Fasta(RefCell::new(Vec::new())),
                accessed,
                modified,
            ),
            seq_file: Fragment::make_virtual_file(
                seq_ino,
                &format!("{}{}", id, SEQ_EXT),
                0o664,
                data_size,
                FileClass::Seq,
                accessed,
                modified,
            ),
        }
    }

    fn rename(&mut self, new_id: &str) {
        self.id = new_id.into();
        self.refresh_virtual_files();
    }

    fn refresh_virtual_files(&mut self) {
        self.fasta_file.name = format!("{}{}", self.id, FASTA_EXT).into();
        self.fasta_file.attrs.size = (self.label_size() + self.data_size()) as u64;

        self.seq_file.name = format!("{}{}", self.id, SEQ_EXT).into();
        self.seq_file.attrs.size = self.data_size() as u64;
    }

    fn label_size(&self) -> usize {
        self.label().len()
    }

    fn data_size(&self) -> usize {
        match &self.data {
            Backing::File(_, start, end) => end - start,
            Backing::Buffer(ref b) | Backing::PureBuffer(ref b) => b.len(),
            Backing::MMap(ref mmap) => mmap.len(),
        }
    }

    fn label(&self) -> String {
        Fragment::make_label(&self.id, &self.name)
    }

    fn data(&self) -> Box<[u8]> {
        match &self.data {
            Backing::File(filename, start, _) => {
                let mut buffer = vec![0u8; self.data_size() as usize];
                let mut f = fs::File::open(filename.as_str())
                    .unwrap_or_else(|_| panic!("Unable to open `{}`", filename));
                f.seek(SeekFrom::Start(*start as u64))
                    .unwrap_or_else(|_| panic!("Unable to seek in `{}`", filename));
                f.read_exact(&mut buffer)
                    .unwrap_or_else(|_| panic!("Unable to read from `{}`", filename));
                buffer.into_boxed_slice()
            }
            Backing::Buffer(ref b) | Backing::PureBuffer(ref b) => b[..].into(),
            Backing::MMap(ref mmap) => mmap[..].into(),
        }
    }

    fn chunk(&self, offset: i64, size: u32) -> Box<[u8]> {
        match &self.data {
            Backing::File(filename, start, _) => {
                let size = size as usize;
                let offset = offset as usize;
                let mut buffer = vec![0u8; size as usize];
                let mut f = fs::File::open(filename.as_str())
                    .unwrap_or_else(|_| panic!("Unable to open `{}`", filename));
                f.seek(SeekFrom::Start(*start as u64 + offset as u64))
                    .unwrap_or_else(|_| panic!("Unable to seek in `{}`", filename));

                f.read_exact(&mut buffer).unwrap_or_else(|e| {
                    panic!(
                        "Unable to read {}:{} from `{}` ({}): {}",
                        start + offset,
                        start + offset + size,
                        filename,
                        f.metadata().unwrap().len(),
                        e
                    )
                });
                buffer.into_boxed_slice()
            }
            Backing::Buffer(ref b) | Backing::PureBuffer(ref b) => {
                let start = offset as usize;
                let end = std::cmp::min(b.len() as i64, offset + size as i64) as usize;
                b[start..end].into()
            }
            Backing::MMap(ref mmap) => {
                let start = offset as usize;
                let end = std::cmp::min(mmap.len() as i64, offset + size as i64) as usize;
                mmap[start..end].into()
            }
        }
    }

    // The same as `chunk`, but skipping new lines.
    fn pure_chunk(&self, offset: i64, size: u32) -> Box<[u8]> {
        let offset = offset as usize;
        let size = size as usize;
        match &self.data {
            Backing::File(filename, start, _) => {
                let mut f = fs::File::open(filename.as_str())
                    .unwrap_or_else(|_| panic!("Unable to open `{}`", filename));
                f.seek(SeekFrom::Start(*start as u64 + offset as u64))
                    .unwrap_or_else(|_| panic!("Unable to seek in `{}`", filename));
                f.bytes()
                    .filter_map(Result::ok)
                    .filter(|&c| c != b'\n')
                    .skip(offset)
                    .take(size)
                    .collect::<Vec<_>>()
                    .into()
            }
            Backing::Buffer(ref b) => b
                .iter()
                .cloned()
                .filter(|&c| c != b'\n')
                .skip(offset)
                .take(size)
                .collect::<Vec<_>>()
                .into(),
            Backing::PureBuffer(ref b) => b[offset..offset + size].to_vec().into(),
            Backing::MMap(ref mmap) => mmap
                .iter()
                .cloned()
                .filter(|&c| c != b'\n')
                .skip(offset as usize)
                .take(size as usize)
                .collect::<Vec<_>>()
                .into(),
        }
    }

    fn extend(&mut self, size: usize) {
        match &mut self.data {
            Backing::Buffer(ref mut b) => b.resize_with(size, Default::default),
            _ => unreachable!(),
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

#[derive(PartialEq, Debug)]
pub enum Cache {
    Mmap, // Store fragments as mmapped-memory
    File, // ...or as filename:start-end pairs
    RAM,  // Buffer all fragments in cache
}
pub struct FustaSettings {
    pub cache: Cache,
    pub concretize_threshold: usize, // How much leeway do we have in memory consumption (in B)
    pub csv_separator: String,
    pub no_overwrite: bool,
}

#[derive(Debug)]
struct PendingAppend {
    data: Vec<u8>,
    attrs: FileAttr,
    fasta_ino: u64,
    seq_ino: u64,
}

/// A Subfragment represents a portion of a fragment (chr:start-end)
#[derive(Debug)]
struct SubFragment {
    /// The parent fragment
    fragment: String,
    /// Starting offset in the parent fragment
    start: isize,
    /// The attributes of this subfragment inode. NOTE: that the subfragment
    /// size is encoded in those.
    attrs: FileAttr,
}
impl SubFragment {
    fn new(fragment: &str, start: isize, attrs: FileAttr) -> SubFragment {
        SubFragment {
            fragment: fragment.to_owned(),
            start,
            attrs,
        }
    }
}
lazy_static! {
    static ref SUBFRAGMENT_RE: Regex = Regex::new(r"(.+):(\d+)-(\d+)").unwrap();
}

pub struct FustaFS {
    fragments: Vec<Fragment>,
    name2fragment: HashMap<String, usize>,
    ino2fragment: HashMap<u64, usize>,

    metadata: fs::Metadata,
    dir_attrs: BTreeMap<u64, FileAttr>,
    files: Vec<Box<dyn VirtualFile + Send>>,
    filename: String,
    settings: FustaSettings,
    current_ino: u64,

    pending_appends: BTreeMap<String, PendingAppend>,

    subfragments: MultiMap<String, u64, SubFragment>, // name -> inode -> SubFragment

    dirty: bool,
}

impl FustaFS {
    pub fn new(settings: FustaSettings, filename: &str) -> Result<FustaFS> {
        let metadata = fs::metadata(filename).context(format!("while opening `{}`", filename))?;
        let mut r = FustaFS {
            fragments: Vec::new(),
            name2fragment: HashMap::new(),
            ino2fragment: HashMap::new(),

            filename: String::new(),
            dir_attrs: btreemap! {
                // Virtual folders
                ROOT_DIR         => FustaFS::make_dir_attrs(ROOT_DIR, 0o775),
                SEQ_DIR          => FustaFS::make_dir_attrs(SEQ_DIR, 0o775),
                FASTA_DIR        => FustaFS::make_dir_attrs(FASTA_DIR, 0o555),
                APPEND_DIR       => FustaFS::make_dir_attrs(APPEND_DIR, 0o775),
                SUBFRAGMENTS_DIR => FustaFS::make_dir_attrs(SUBFRAGMENTS_DIR, 0o555),
            },
            files: vec![
                Box::new(BufferFile {
                    name: INFO_FILE_NAME.into(),
                    ino: INFO_FILE,
                    attrs: FustaFS::make_file_attrs(INFO_FILE, 0o444),
                    class: FileClass::Text,
                    _data: Vec::new(),
                }),
                Box::new(BufferFile {
                    name: INFO_CSV_FILE_NAME.into(),
                    ino: INFO_CSV_FILE,
                    attrs: FustaFS::make_file_attrs(INFO_CSV_FILE, 0o444),
                    class: FileClass::Text,
                    _data: Vec::new(),
                }),
                Box::new(BufferFile {
                    name: LABELS_FILE_NAME.into(),
                    ino: LABELS_FILE,
                    attrs: FustaFS::make_file_attrs(LABELS_FILE, 0o444),
                    class: FileClass::Text,
                    _data: Vec::new(),
                }),
            ],
            metadata,
            settings,
            current_ino: FIRST_INO,
            pending_appends: Default::default(),
            subfragments: Default::default(),
            dirty: false,
        };

        r.read_fasta(filename)
            .context(format!("while parsing {}", filename))?;
        Ok(r)
    }

    fn make_dir_attrs(ino: u64, perms: u16) -> FileAttr {
        FileAttr {
            ino,
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
            blksize: 0,
        }
    }

    fn make_file_attrs(ino: u64, perms: u16) -> FileAttr {
        FileAttr {
            ino,
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
            blksize: 512,
        }
    }

    fn make_file_attrs_with_size(ino: u64, perms: u16, size: u64) -> FileAttr {
        let mut attrs = Self::make_file_attrs(ino, perms);
        attrs.size = size;
        attrs
    }

    fn new_ino(&mut self) -> u64 {
        let r = self.current_ino;
        self.current_ino += 1;
        r
    }

    fn get_file(&mut self, ino: u64) -> Option<&mut Box<dyn VirtualFile + Send>> {
        self.files.iter_mut().find(|f| f.ino() == ino)
    }

    fn read_fasta(&mut self, filename: &str) -> Result<()> {
        info!("Reading {}...", filename);
        let fasta_file =
            fs::File::open(filename).context(format!("failed to open file `{}`", filename))?;
        let fragments =
            FastaReader::new(fasta_file, self.settings.cache == Cache::RAM).collect::<Vec<_>>();
        let mut keys = fragments.iter().map(|f| &f.id).collect::<Vec<_>>();
        keys.sort();
        keys.dedup();
        if keys.len() != fragments.len() {
            anyhow::bail!("`{}` contains duplicated keys", filename)
        }

        let file = fs::File::open(filename).context(format!("Failed to open `{}`", filename))?;
        self.filename = filename.to_owned();

        self.fragments = fragments
            .into_iter()
            .map(|fragment| {
                if fragment.id.chars().any(|c| FORBIDDEN_CHARS.contains(&c)) {
                    Err(anyhow::anyhow!(format!(
                        "Fragment ID `{}` contains a forbidden character",
                        fragment.id
                    )))
                } else {
                    Ok(Fragment::new(
                        &fragment.id,
                        &fragment.name,
                        match self.settings.cache {
                            Cache::Mmap => Backing::MMap(unsafe {
                                memmap2::MmapOptions::new()
                                    .offset(fragment.pos.0 as u64)
                                    .len(fragment.len)
                                    .map(&file)
                                    .unwrap()
                            }),
                            Cache::File => {
                                Backing::File(filename.into(), fragment.pos.0, fragment.pos.1)
                            }
                            Cache::RAM => Backing::PureBuffer(fragment.seq.unwrap()),
                        },
                        self.new_ino(),
                        self.new_ino(),
                        self.metadata.accessed().unwrap(),
                        self.metadata.modified().unwrap(),
                    ))
                }
            })
            .collect::<Result<Vec<_>>>()?;
        self.refresh_metadata(true);
        info!("Done.");
        Ok(())
    }

    fn concretize(&mut self, force: bool) {
        if !self.dirty {
            debug!("CONCRETIZE: nothing to do");
            return;
        }

        let in_memory = self.fragments.iter().fold(0, |ax, f| {
            ax + match &f.data {
                Backing::Buffer(b) => b.len(),
                _ => 0,
            }
        });

        // We only concretize if the call is not forced and
        // 1. the allowed cache is not yet used
        // or
        // 2. the user wants to cache everything anyway
        if !force
            && (in_memory < self.settings.concretize_threshold
                || !(self.settings.cache == Cache::RAM))
        {
            trace!(
                "Using {:.2}MB out of {}; skipping concretization",
                in_memory / (1024 * 1024),
                self.settings.concretize_threshold / (1024 * 1024)
            );
            return;
        }

        trace!("========== CONCRETIZING ========");
        notify(format!("Updating {}", &self.filename));
        let mut index = 0;
        let mut last_start;
        trace!("Writing fragments");
        let tmp_filename = format!("{}#fusta#", &self.filename);
        {
            // Scope to ensure the tmp file is correctly closed
            let mut tmp_file = fs::File::create(&tmp_filename)
                .unwrap_or_else(|_| panic!("Unable to create `{}`", tmp_filename));
            for fragment in self.fragments.iter_mut() {
                trace!("Writing {}", fragment.id);
                tmp_file
                    .write_all(fragment.label().as_bytes())
                    .unwrap_or_else(|_| panic!("Unable to write to `{}`", tmp_filename));
                index += fragment.label().as_bytes().len();
                last_start = index;
                tmp_file
                    .write_all(&fragment.data())
                    .unwrap_or_else(|_| panic!("Unable to write to `{}`", tmp_filename));
                index += fragment.data().len();
                if let Some(c) = fragment.data().last() {
                    if *c != b'\n' {
                        tmp_file
                            .write_all(&[b'\n'])
                            .unwrap_or_else(|_| panic!("Unable to write to `{}`", tmp_filename));
                        index += 1;
                    }
                }

                fragment.data = Backing::File(self.filename.clone().into(), last_start, index);
                fragment.refresh_virtual_files();
            }
        }
        trace!("Renaming {} to {}", tmp_filename, &self.filename);
        fs::rename(&tmp_filename, &self.filename).unwrap_or_else(|_| {
            panic!(
                "Unable to rename `{}` to `{}`",
                &tmp_filename, &self.filename
            )
        });
        trace!("========== DONE ========");
        notify(format!("{} has been updated", &self.filename));
        self.dirty = false;
    }

    fn fragment_from_id(&self, id: &str) -> Option<&Fragment> {
        self.name2fragment.get(id).map(|&i| &self.fragments[i])
    }

    fn fragment_from_ino(&self, ino: u64) -> Option<&Fragment> {
        self.ino2fragment
            .get(&ino)
            .and_then(|&i| self.fragments.get(i))
    }

    fn mut_fragment_from_ino(&mut self, ino: u64) -> Option<&mut Fragment> {
        if let Some(i) = self.ino2fragment.get(&ino) {
            self.fragments.get_mut(*i)
        } else {
            None
        }
    }

    fn fragment_from_fasta_filename(&self, name: &str) -> Option<&Fragment> {
        name.strip_suffix(FASTA_EXT)
            .and_then(|id| self.name2fragment.get(id))
            .and_then(|i| self.fragments.get(*i))
    }

    fn fragment_from_seq_filename(&self, name: &str) -> Option<&Fragment> {
        name.strip_suffix(SEQ_EXT)
            .and_then(|id| self.name2fragment.get(id))
            .and_then(|i| self.fragments.get(*i))
    }

    fn mut_fragment_from_seq_filename(&mut self, name: &str) -> Option<&mut Fragment> {
        let id = name.strip_suffix(SEQ_EXT).unwrap();
        if let Some(i) = self.name2fragment.get(id) {
            Some(&mut self.fragments[*i])
        } else {
            None
        }
    }

    fn mut_fragment_from_fasta_filename(&mut self, name: &str) -> Option<&mut Fragment> {
        let id = name.strip_suffix(FASTA_EXT).unwrap();
        if let Some(i) = self.name2fragment.get(id) {
            Some(&mut self.fragments[*i])
        } else {
            None
        }
    }

    fn subfragment_from_ino(&self, ino: u64) -> Option<&SubFragment> {
        self.subfragments.get_alt(&ino)
    }

    fn make_info_buffer(&mut self) {
        use ascii_table::*;
        use num_format::*;

        let mut table = AsciiTable::default();
        table.column(0).set_header("ID");
        table.column(1).set_header("Infos");
        table.column(2).set_header("Length (bp)");

        trace!("Making INFO BUFFER");
        let header = format!("{} - {} sequences", &self.filename, self.fragments.len());
        let infos = self
            .fragments
            .iter()
            .map(|f| {
                vec![
                    f.id.to_owned(),
                    f.name.as_ref().unwrap_or(&"".to_string()).into(),
                    f.data_size().to_formatted_string(&Locale::en).into(),
                ]
            })
            .collect::<Vec<_>>();
        let content = format!(
            "{}\n{}\n{}",
            header,
            "=".repeat(header.len()),
            &table.format(&infos)
        );
        let size = content.as_bytes().len() as u64;
        if let Some(x) = self.get_file(INFO_FILE) {
            x.set_data(content.as_bytes());
            x.mut_attrs().size = size;
        }
    }

    fn make_info_csv_buffer(&mut self) {
        trace!("Making INFO_CSV BUFFER");
        let header = format!(
            "id{}name{}length",
            self.settings.csv_separator, self.settings.csv_separator
        );
        let infos = self
            .fragments
            .iter()
            .map(|f| {
                format!(
                    "{}{}\"{}\"{}{}",
                    f.id.to_owned(),
                    self.settings.csv_separator,
                    f.name.as_ref().unwrap_or(&"".to_string()),
                    self.settings.csv_separator,
                    f.data_size()
                )
            })
            .collect::<Vec<_>>();
        let content = format!("{}\n{}", header, infos.join("\n"));
        let size = content.as_bytes().len() as u64;
        if let Some(x) = self.get_file(INFO_CSV_FILE) {
            x.set_data(content.as_bytes());
            x.mut_attrs().size = size;
        }
    }

    fn make_labels_buffer(&mut self) {
        trace!("Making LABELS BUFFER");
        let content = self
            .fragments
            .iter()
            .map(Fragment::label)
            .collect::<Vec<_>>()
            .join("");
        let size = content.as_bytes().len() as u64;
        if let Some(x) = self.get_file(LABELS_FILE) {
            x.set_data(content.as_bytes());
            x.mut_attrs().size = size;
        }
    }

    fn refresh_metadata(&mut self, force: bool) {
        if self.dirty || force {
            debug!("Refreshing metadata...");
            self.make_info_buffer();
            self.make_info_csv_buffer();
            self.make_labels_buffer();
            self.update_indices();
            debug!("Done.")
        }
    }

    fn update_indices(&mut self) {
        self.name2fragment = self
            .fragments
            .iter()
            .enumerate()
            .map(|(i, f)| (f.id.clone().into(), i))
            .collect::<HashMap<_, _>>();

        self.ino2fragment = self
            .fragments
            .iter()
            .enumerate()
            .flat_map(|(i, f)| {
                std::iter::IntoIterator::into_iter([(f.fasta_file.ino, i), (f.seq_file.ino, i)])
            })
            .collect::<HashMap<_, _>>();
    }

    fn is_fasta_file(&self, ino: u64) -> bool {
        self.ino2fragment
            .get(&ino)
            .and_then(|i| self.fragments.get(*i))
            .map(|f| f.fasta_file.ino == ino)
            .unwrap_or(false)
    }

    fn is_seq_file(&self, ino: u64) -> bool {
        self.ino2fragment
            .get(&ino)
            .and_then(|i| self.fragments.get(*i))
            .map(|f| f.seq_file.ino == ino)
            .unwrap_or(false)
    }

    fn is_append_file(&self, ino: u64) -> bool {
        self.pending_appends.iter().any(|p| p.1.attrs.ino == ino)
    }

    fn is_writeable(&self, ino: u64) -> bool {
        self.is_append_file(ino) || self.is_seq_file(ino)
    }

    fn create_subfragment(&mut self, name: &str) -> Result<FileAttr, String> {
        fn clear_coordinates(start: isize, end: isize) -> (isize, isize) {
            let start = if start < 0 {
                warn!("Invalid start position {}; using 1 instead", start);
                0
            } else {
                start
            };
            let end = if end < start {
                warn!(
                    "{}-{}: {} < {}; using {} instead",
                    start,
                    end,
                    end,
                    start,
                    std::cmp::max(start, end)
                );
                std::cmp::max(start, end)
            } else {
                end
            };
            (start, end)
        }

        let error_message = format!("`{}` is not a valid subfragment scheme", name);
        if name.contains(':') && name.contains('-') {
            let caps = SUBFRAGMENT_RE.captures(name).ok_or_else(|| {
                format!("{}: it should be of the form ID:START-END", error_message)
            })?;
            if caps.len() == 4 {
                let ino = self.new_ino();
                let fragment_id = self
                    .fragment_from_id(&caps[1])
                    .ok_or_else(|| format!("`{}` is not a fragment", &caps[1]))?
                    .id
                    .clone();
                let start = str::parse::<isize>(&caps[2])
                    .map_err(|_| format!("{}: `{}` is not an integer", &error_message, &caps[1]))?
                    - 1;
                let end = str::parse::<isize>(&caps[3])
                    .map_err(|_| format!("{}: `{}` is not an integer", &error_message, &caps[2]))?;
                let (start, end) = clear_coordinates(start, end);
                let key = format!("{}:{}-{}", fragment_id, start, end);

                self.subfragments
                    .get(&key)
                    .map(|sf| sf.attrs)
                    .or_else(|| {
                        let size = (end - start) as u64;
                        let attrs = FustaFS::make_file_attrs_with_size(ino, 0o444, size);
                        let sf = SubFragment::new(&fragment_id, start, attrs);
                        self.subfragments.insert(key, ino, sf);
                        Some(attrs)
                    })
                    .ok_or_else(|| unimplemented!())
            } else {
                Err(error_message)
            }
        } else {
            let ino = self.new_ino();
            let fragment = self
                .fragment_from_id(name)
                .ok_or_else(|| format!("`{}` is not a fragment", name))?;
            let fragment_id = fragment.id.clone();
            let fragment_len = fragment.data.len();
            let attrs =
                FustaFS::make_file_attrs_with_size(ino, 0o444, fragment_len.try_into().unwrap());
            let sf = SubFragment::new(&fragment_id, 0, attrs);
            self.subfragments.insert(name.to_owned(), ino, sf);
            Ok(attrs)
        }
    }
}

impl Drop for FustaFS {
    fn drop(&mut self) {
        self.concretize(true);
    }
}
impl Filesystem for FustaFS {
    fn lookup(&mut self, _req: &Request, parent: u64, name: &OsStr, reply: ReplyEntry) {
        let name = name.to_str().unwrap();
        match parent {
            ROOT_DIR => match name {
                "fasta" => {
                    reply.entry(&TTL, &self.dir_attrs[&FASTA_DIR], 0);
                }
                "seqs" => {
                    reply.entry(&TTL, &self.dir_attrs[&SEQ_DIR], 0);
                }
                "append" => {
                    reply.entry(&TTL, &self.dir_attrs[&APPEND_DIR], 0);
                }
                "get" => {
                    reply.entry(&TTL, &self.dir_attrs[&SUBFRAGMENTS_DIR], 0);
                }
                INFO_FILE_NAME => {
                    reply.entry(&TTL, self.get_file(INFO_FILE).unwrap().attrs(), 0);
                }
                INFO_CSV_FILE_NAME => {
                    reply.entry(&TTL, self.get_file(INFO_CSV_FILE).unwrap().attrs(), 0);
                }
                LABELS_FILE_NAME => {
                    reply.entry(&TTL, self.get_file(LABELS_FILE).unwrap().attrs(), 0);
                }
                _ => {
                    reply.error(ENOENT);
                }
            },
            SEQ_DIR | FASTA_DIR => {
                let file = (if parent == SEQ_DIR {
                    self.fragment_from_seq_filename(name)
                } else {
                    self.fragment_from_fasta_filename(name)
                })
                .and_then(|f| f.file_from_filename(name));

                if let Some(file) = file {
                    reply.entry(&TTL, file.attrs(), 0);
                } else {
                    reply.error(ENOENT);
                }
            }
            SUBFRAGMENTS_DIR => {
                let sf = self.create_subfragment(name);
                match sf {
                    Ok(attrs) => {
                        reply.entry(&TTL, &attrs, 0);
                    }
                    Err(e) => {
                        warn!("{}", &e);
                        reply.error(ENOENT);
                    }
                }
            }
            _ => {
                warn!("LOOKUP: parent {} does not exist", parent);
                reply.error(ENOENT);
            }
        }
    }

    fn getattr(&mut self, _req: &Request, ino: u64, reply: ReplyAttr) {
        match ino {
            ROOT_DIR | SEQ_DIR | APPEND_DIR | FASTA_DIR | SUBFRAGMENTS_DIR => {
                reply.attr(&TTL, self.dir_attrs.get(&ino).unwrap())
            }
            INFO_FILE => reply.attr(&TTL, self.get_file(INFO_FILE).unwrap().attrs()),
            INFO_CSV_FILE => reply.attr(&TTL, self.get_file(INFO_CSV_FILE).unwrap().attrs()),
            LABELS_FILE => reply.attr(&TTL, self.get_file(LABELS_FILE).unwrap().attrs()),
            ino if self.subfragment_from_ino(ino).is_some() => {
                reply.attr(&TTL, &self.subfragment_from_ino(ino).unwrap().attrs);
            }
            _ => {
                if let Some(file) = self
                    .fragment_from_ino(ino)
                    .and_then(|f| f.file_from_ino(ino))
                {
                    reply.attr(&TTL, file.attrs())
                } else {
                    warn!("GETATTR: ino `{}` does not exist", ino);
                    reply.error(ENOENT)
                }
            }
        }
    }

    fn read(
        &mut self,
        _req: &Request,
        ino: u64,
        _fh: u64,
        offset: i64,
        size: u32,
        _flags: i32,
        _lock_owner: Option<u64>,
        reply: ReplyData,
    ) {
        debug!("READING {}", ino);
        match ino {
            INFO_FILE => {
                let data = self.get_file(INFO_FILE).unwrap().data();
                let start = offset as usize;
                let end = std::cmp::min(start + size as usize, data.len());
                reply.data(&data[start..end]);
            }
            INFO_CSV_FILE => {
                let data = self.get_file(INFO_CSV_FILE).unwrap().data();
                let start = offset as usize;
                let end = std::cmp::min(start + size as usize, data.len());
                reply.data(&data[start..end]);
            }
            LABELS_FILE => {
                let data = self.get_file(LABELS_FILE).unwrap().data();
                let start = offset as usize;
                let end = std::cmp::min(start + size as usize, data.len());
                reply.data(&data[start..end]);
            }
            ino if self.ino2fragment.contains_key(&ino) => {
                let fragment = match self.mut_fragment_from_ino(ino) {
                    Some(f) => f,
                    _ => {
                        error!("No fragment linked to ino {}", ino);
                        return;
                    }
                };
                match fragment
                    .file_from_ino(ino)
                    .expect("No file linked to this fragment")
                    .class()
                {
                    FileClass::Fasta(_) => {
                        let label_size = fragment.label_size() as i64;
                        if offset > label_size as i64 {
                            reply.data(&fragment.chunk(offset - label_size, size))
                        } else {
                            let end = offset as usize + size as usize;
                            let label = fragment.label();
                            let data_chunk =
                                fragment.chunk(0, std::cmp::min(fragment.data_size() as u32, size));
                            match fragment
                                .mut_file_from_ino(ino)
                                .expect("No file linked to this fragment")
                                .class()
                            {
                                FileClass::Fasta(header_buffer) => {
                                    if end > header_buffer.borrow().len() {
                                        header_buffer.replace(
                                            label
                                                .as_bytes()
                                                .iter()
                                                .chain(data_chunk.iter())
                                                .cloned()
                                                .take(end as usize)
                                                .collect::<Vec<_>>(),
                                        );
                                    }
                                    reply.data(
                                        &header_buffer.borrow()[offset as usize
                                            ..std::cmp::min(end, header_buffer.borrow().len())],
                                    );
                                }
                                _ => panic!("WTF"),
                            }
                        }
                    }
                    FileClass::Seq => reply.data(&fragment.chunk(offset, size)),
                    FileClass::Text => unimplemented!(), // A fragment can never refer to a text file
                }
            }
            ino if self.subfragment_from_ino(ino).is_some() => {
                let subfragment = self.subfragment_from_ino(ino).unwrap();
                match self.fragment_from_id(&subfragment.fragment) {
                    Some(fragment) => {
                        reply.data(&fragment.pure_chunk(offset + subfragment.start as i64, size))
                    }
                    _ => {
                        error!("No fragment linked to ino {}", ino);
                        reply.error(ENOENT);
                    }
                };
            }
            _ => {
                warn!("READ: {} is not a file", ino);
                reply.error(ENOENT);
            }
        }
    }

    fn readdir(
        &mut self,
        _req: &Request,
        ino: u64,
        _fh: u64,
        offset: i64,
        mut reply: ReplyDirectory,
    ) {
        match ino {
            ROOT_DIR => {
                let entries = btreemap! {
                    ROOT_DIR         => (FileType::Directory, "."),
                    0                => (FileType::Directory, ".."), // TODO
                    FASTA_DIR        => (FileType::Directory, "fasta"),
                    SEQ_DIR          => (FileType::Directory, "seqs"),
                    APPEND_DIR       => (FileType::Directory, "append"),
                    SUBFRAGMENTS_DIR => (FileType::Directory, "get"),
                    INFO_FILE        => (FileType::RegularFile, INFO_FILE_NAME),
                    INFO_CSV_FILE    => (FileType::RegularFile, INFO_CSV_FILE_NAME),
                    LABELS_FILE      => (FileType::RegularFile, LABELS_FILE_NAME),
                };
                for (o, (ino, entry)) in entries.iter().enumerate().skip(offset as usize) {
                    let _ = reply.add(*ino, o as i64 + 1, entry.0, &entry.1);
                }
                reply.ok();
            }
            FASTA_DIR | SEQ_DIR => {
                for (i, file) in vec![
                    (ino, FileType::Directory, &".".into()),
                    (ROOT_DIR, FileType::Directory, &"..".into()),
                ]
                .into_iter()
                .chain(self.fragments.iter().map(|f| {
                    let (ino, name) = if ino == SEQ_DIR {
                        (f.seq_file.ino, &f.seq_file.name)
                    } else {
                        (f.fasta_file.ino, &f.fasta_file.name)
                    };
                    (ino, FileType::RegularFile, name)
                }))
                .enumerate()
                .skip(offset.try_into().unwrap())
                {
                    if reply.add(file.0, (i + 1).try_into().unwrap(), file.1, file.2.as_str()) {
                        break;
                    }
                }
                reply.ok();
            }
            APPEND_DIR => {
                let entries = btreemap! {
                    APPEND_DIR => (FileType::Directory, ".".to_owned()),
                    ROOT_DIR   => (FileType::Directory, "..".to_owned()),
                };
                for (o, (ino, entry)) in entries.iter().enumerate().skip(offset as usize) {
                    let _ = reply.add(*ino, o as i64 + 1, entry.0, &entry.1);
                }
                reply.ok();
            }
            SUBFRAGMENTS_DIR => {
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
            }
            SEQ_DIR | FASTA_DIR => {
                let name = name.to_str().unwrap();
                let file = (if parent == SEQ_DIR {
                    self.fragment_from_seq_filename(name)
                } else {
                    self.fragment_from_fasta_filename(name)
                })
                .and_then(|f| f.file_from_filename(name));

                if file.is_some() {
                    let length_before = self.fragments.len();
                    self.fragments
                        .retain(|f| f.fasta_file.name != name && f.seq_file.name != name);
                    let length_after = self.fragments.len();
                    // Only mark as dirty if we effectively removed something
                    if length_after != length_before {
                        self.dirty = true;
                    }
                    self.refresh_metadata(false);
                    self.concretize(false);
                    reply.ok();
                } else {
                    warn!("UNLINK: unknown file: `{:?}`", name);
                    reply.error(ENOENT);
                }
            }
            APPEND_DIR => {
                warn!("UNLINK: unauthorized in append virtual dir");
                reply.error(EACCES);
            }
            SUBFRAGMENTS_DIR => {
                reply.error(ENOENT);
            }
            _ => {
                warn!("UNLINK: parent {} does not exist", parent);
                reply.error(ENOENT);
            }
        }
    }

    fn mknod(
        &mut self,
        _req: &Request,
        parent: u64,
        name: &OsStr,
        _mode: u32,
        _umask: u32,
        _rdev: u32,
        reply: ReplyEntry,
    ) {
        match parent {
            ROOT_DIR | SEQ_DIR | FASTA_DIR | SUBFRAGMENTS_DIR => {
                warn!("MKNOD: writing in {} is forbidden", parent);
                reply.error(EACCES);
            }
            APPEND_DIR => {
                let name = name.to_str().unwrap();
                let basename = std::path::Path::new(name)
                    .file_stem()
                    .unwrap()
                    .to_str()
                    .unwrap();
                // From man: if pathname already exists [...], this call fails with an EEXIST error.
                if self.fragments.iter().any(|f| f.id == basename) && self.settings.no_overwrite {
                    notify(format!("Cannot create `{:?}`, already exists", name));
                    error!("Cannot create `{:?}`, already exists", name);
                    reply.error(EEXIST);
                    return;
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
                    blksize: 512,
                };
                let pending = PendingAppend {
                    data: Vec::new(),
                    attrs: attrs,
                    seq_ino: self.new_ino(),
                    fasta_ino: self.new_ino(),
                };
                reply.entry(&TTL, &pending.attrs, 0);
                self.pending_appends.insert(basename.to_string(), pending);
            }
            _ => {
                warn!("MKNOD: parent {} does not exist", parent);
                reply.error(ENOENT);
            }
        }
    }

    fn write(
        &mut self,
        _req: &Request,
        ino: u64,
        _fh: u64,
        offset: i64,
        data: &[u8],
        _write_flags: u32,
        _flags: i32,
        _lock_owner: Option<u64>,
        reply: ReplyWrite,
    ) {
        if !self.is_writeable(ino) {
            error!("{} is not writeable", ino);
            reply.error(EACCES);
        } else {
            // Check that we only write valid FASTA content
            if !data.iter().all(|&c| is_fasta_char(c)) {
                warn!("Invalid FASTA sequence found when writing");
                notify(format!("Invalid FASTA sequence found when writing"));
                reply.error(EIO);
                return;
            }

            // We write to an existing fragment
            if self.fragment_from_ino(ino).is_some() {
                let fragment = self
                    .mut_fragment_from_ino(ino)
                    .expect("Something went very wrong");
                // As soon as there's a write, we have to switch this fragment to a buffer-backed storage
                if !matches!(fragment.data, Backing::Buffer(_)) {
                    fragment.data = Backing::Buffer(fragment.data().to_vec());
                }

                // Ensure that the backing buffer is big enough
                let max_size = offset as usize + data.len();
                if max_size > fragment.data_size() {
                    fragment.extend(max_size)
                }

                // Then finally write the data
                if let Backing::Buffer(b) = &mut fragment.data {
                    let start = offset as usize;
                    let end = start + data.len();
                    b.splice(start..end, data.iter().cloned());
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
                pending_fragment
                    .data
                    .splice(start..end, data.iter().cloned());
                reply.written(data.len() as u32);
                trace!("\tWriting to {}", name);
            } else {
                reply.error(ENOENT);
            }
        }
    }

    fn setattr(
        &mut self,
        _req: &Request<'_>,
        ino: u64,
        mode: Option<u32>,
        uid: Option<u32>,
        gid: Option<u32>,
        size: Option<u64>,
        atime: Option<TimeOrNow>,
        mtime: Option<TimeOrNow>,
        _ctime: Option<SystemTime>,
        _fh: Option<u64>,
        crtime: Option<SystemTime>,
        chgtime: Option<SystemTime>,
        bkuptime: Option<SystemTime>,
        flags: Option<u32>,
        reply: ReplyAttr,
    ) {
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
            ROOT_DIR | SEQ_DIR | FASTA_DIR => reply.error(EACCES),
            INFO_FILE | INFO_CSV_FILE | LABELS_FILE => reply.error(EACCES),
            _ => {
                if self.fragment_from_ino(ino).is_some() {
                    if !self.is_writeable(ino) {
                        reply.error(EACCES);
                    } else {
                        if let Some(file) = self
                            .mut_fragment_from_ino(ino)
                            .and_then(|f| f.mut_file_from_ino(ino))
                        {
                            if let Some(uid) = uid {
                                file.mut_attrs().uid = uid
                            }
                            if let Some(gid) = gid {
                                file.mut_attrs().gid = gid
                            }
                            // if let Some(atime) = atime {
                            //     file.mut_attrs().atime = atime
                            // }
                            // if let Some(mtime) = mtime {
                            //     file.mut_attrs().mtime = mtime
                            // }
                            if let Some(chgtime) = chgtime {
                                file.mut_attrs().mtime = chgtime
                            }
                            if let Some(crtime) = crtime {
                                file.mut_attrs().crtime = crtime
                            }
                            // macOS only
                            if let Some(flags) = flags {
                                file.mut_attrs().flags = flags
                            }
                            if let Some(mode) = mode {
                                file.mut_attrs().perm = mode as u16
                            }
                        }
                        if let Some(size) = size.map(|s| s as usize) {
                            if size == 0 {
                                // Clear the file, called by the truncate syscall
                                if let Some(fragment) = self.mut_fragment_from_ino(ino) {
                                    fragment.data = Backing::Buffer(Vec::new());
                                    fragment.refresh_virtual_files();
                                }
                            } else if size
                                != self
                                    .fragment_from_ino(ino)
                                    .map(Fragment::data_size)
                                    .unwrap()
                            {
                                // Redim the file
                                let fragment = self.mut_fragment_from_ino(ino).unwrap();
                                if !matches!(fragment.data, Backing::Buffer(_)) {
                                    fragment.data = Backing::Buffer(fragment.data().to_vec());
                                }
                                fragment.extend(size);
                                fragment.refresh_virtual_files();
                                self.dirty = true;
                            }
                            self.mut_fragment_from_ino(ino)
                                .and_then(|f| f.mut_file_from_ino(ino))
                                .unwrap()
                                .mut_attrs()
                                .size = size as u64;
                        }
                        reply.attr(
                            &TTL,
                            self.fragment_from_ino(ino)
                                .and_then(|f| f.file_from_ino(ino))
                                .unwrap()
                                .attrs(),
                        );
                    }
                } else if let Some((name, pending_fragment)) = self
                    .pending_appends
                    .iter_mut()
                    .find(|(_, p)| p.attrs.ino == ino)
                {
                    if let Some(size) = size {
                        trace!("\tResizing {} @{}", name, size);
                        pending_fragment
                            .data
                            .resize_with(size as usize, Default::default);
                        reply.attr(&TTL, &pending_fragment.attrs);
                    }
                } else {
                    warn!("\t{:?} does not exist", ino);
                    reply.error(ENOENT);
                }
            }
        }
    }

    fn destroy(&mut self) {
        info!("Closing FUSTA");
        self.concretize(false)
    }

    fn rename(
        &mut self,
        _req: &Request,
        parent: u64,
        name: &OsStr,
        newparent: u64,
        newname: &OsStr,
        _flags: u32,
        reply: ReplyEmpty,
    ) {
        match parent {
            ROOT_DIR | APPEND_DIR | SUBFRAGMENTS_DIR => {
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
                    if newname_path
                        .extension()
                        .map(|ext| ext == DOTLESS_SEQ_EXT)
                        .unwrap_or(false)
                    {
                        debug!("Using {:?} for {}", newname_path.file_stem(), new_id);
                        new_id = newname_path
                            .file_stem()
                            .unwrap()
                            .to_str()
                            .unwrap()
                            .to_string();
                    }
                    // Shortcut if we cannot overwrite existing fragments
                    let replaced_fragment = self.fragment_from_id(&new_id);
                    if replaced_fragment.is_some() && self.settings.no_overwrite {
                        error!("Cannot rename {:?} to {}: already existing.", name, new_id);
                        reply.error(EACCES);
                    } else {
                        if replaced_fragment.is_some() {
                            warn!("Replacing {}", new_id);
                            self.fragments.retain(|f| f.id != new_id);
                        }
                        if let Some(ref mut fragment) = if parent == SEQ_DIR {
                            self.mut_fragment_from_seq_filename(name.to_str().unwrap())
                        } else {
                            self.mut_fragment_from_fasta_filename(name.to_str().unwrap())
                        } {
                            fragment.rename(&new_id);
                            info!("Renaming {:?} -> {:?}", name, newname);
                            self.dirty = true;
                            self.concretize(false);
                            self.refresh_metadata(false);
                            reply.ok();
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
        self.concretize(false);
        self.refresh_metadata(false);
        reply.ok();
    }

    fn fsyncdir(
        &mut self,
        _req: &Request,
        _ino: u64,
        _fh: u64,
        _datasync: bool,
        reply: ReplyEmpty,
    ) {
        trace!("FSYNCDIR");
        self.concretize(false);
        self.refresh_metadata(false);
        reply.ok();
    }

    fn flush(&mut self, _req: &Request, _ino: u64, _fh: u64, _lock_owner: u64, reply: ReplyEmpty) {
        trace!("FLUSH");
        // self.concretize();
        reply.ok();
    }

    fn release(
        &mut self,
        _req: &Request,
        ino: u64,
        _fh: u64,
        _flags: i32,
        _lock_owner: Option<u64>,
        _flush: bool,
        reply: ReplyEmpty,
    ) {
        debug!("RELEASE {}", ino);
        if self.is_writeable(ino) {
            for pending in self.pending_appends.iter() {
                // FS is dirty at the first pending fragment
                self.dirty = true;

                if pending.1.attrs.ino == ino {
                    trace!("RELEASE: {}", pending.0);
                    trace!("Dumping...");
                    // Dump the pending in its own file and parse it
                    let mut tmpfile =
                        tempfile::tempfile().expect("Unable to create a temporary file");
                    tmpfile
                        .write_all(&pending.1.data)
                        .expect("Unable to write to temporary file");
                    trace!("Parsing...");
                    tmpfile
                        .seek(SeekFrom::Start(0))
                        .expect("Unable to seek in temporary file");
                    let new_fragments = FastaReader::new(&tmpfile, true).collect::<Vec<_>>();
                    let new_keys = new_fragments.iter().map(|f| &f.id).collect::<Vec<_>>();
                    let old_keys = self
                        .fragments
                        .iter()
                        .map(|f| &f.id)
                        .cloned()
                        .collect::<Vec<_>>();

                    // Remove the fragments sharing an existing key if overwrite is allowed
                    let no_overwrite = self.settings.no_overwrite;
                    if !no_overwrite {
                        self.fragments.retain(|f| !new_keys.contains(&&f.id))
                    }
                    self.fragments
                        .extend(new_fragments.into_iter().filter_map(|new_fragment| {
                            if old_keys.contains(&new_fragment.id) && no_overwrite {
                                error!("Skipping `{}`, already existing", &new_fragment.id);
                                None
                            } else {
                                Some(Fragment::new(
                                    &new_fragment.id,
                                    &new_fragment.name,
                                    Backing::PureBuffer(new_fragment.seq.unwrap()),
                                    pending.1.seq_ino,
                                    pending.1.fasta_ino,
                                    pending.1.attrs.atime,
                                    pending.1.attrs.mtime,
                                ))
                            }
                        }));
                }
            }
            self.concretize(false);
            self.refresh_metadata(false);
        } else {
            debug!("Not a writeable file; ignoring")
        }
        reply.ok();
    }
}
