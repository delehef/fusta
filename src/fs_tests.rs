use super::*;
use std::io::Write;
use tempfile::NamedTempFile;

fn fasta(content: &str) -> NamedTempFile {
    let mut f = NamedTempFile::new().unwrap();
    f.write_all(content.as_bytes()).unwrap();
    f.flush().unwrap();
    f
}

fn ram_settings() -> FustaSettings {
    FustaSettings {
        cache: Cache::RAM,
        concretize_threshold: 500 * 1024 * 1024,
        csv_separator: ",".to_string(),
        no_overwrite: false,
    }
}

fn file_settings() -> FustaSettings {
    FustaSettings {
        cache: Cache::File,
        ..ram_settings()
    }
}

fn mmap_settings() -> FustaSettings {
    FustaSettings {
        cache: Cache::Mmap,
        ..ram_settings()
    }
}

#[test]
fn new_with_single_fragment() {
    let f = fasta(">chr1\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(fs.fragments.len(), 1);
}

#[test]
fn new_with_multiple_fragments() {
    let f = fasta(">chr1\nACGT\n>chr2\nGGGG\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(fs.fragments.len(), 2);
}

#[test]
fn new_with_empty_fasta_succeeds() {
    let f = fasta("");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(fs.fragments.len(), 0);
}

#[test]
fn new_nonexistent_file_errors() {
    assert!(FustaFS::new(ram_settings(), "/nonexistent/path/file.fa").is_err());
}

#[test]
fn new_duplicate_ids_errors() {
    let f = fasta(">chr1\nACGT\n>chr1\nGGGG\n");
    assert!(FustaFS::new(ram_settings(), f.path().to_str().unwrap()).is_err());
}

#[test]
fn new_null_char_in_id_errors() {
    let f = fasta(">\0bad\nACGT\n");
    assert!(FustaFS::new(ram_settings(), f.path().to_str().unwrap()).is_err());
}

#[cfg(unix)]
#[test]
fn new_backslash_in_id_errors() {
    let f = fasta(">bad\\id\nACGT\n");
    assert!(FustaFS::new(ram_settings(), f.path().to_str().unwrap()).is_err());
}

#[test]
fn lookup_by_id_found_and_not_found() {
    let f = fasta(">chr1\nAAAA\n>chr2\nCCCC\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert!(fs.fragment_from_id("chr1").is_some());
    assert!(fs.fragment_from_id("chr2").is_some());
    assert!(fs.fragment_from_id("chr99").is_none());
}

#[test]
fn lookup_by_fasta_ino() {
    let f = fasta(">chr1\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    let ino = fs.fragments[0].fasta_file.ino;
    assert!(fs.fragment_from_ino(ino).is_some());
}

#[test]
fn lookup_by_seq_ino() {
    let f = fasta(">chr1\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    let ino = fs.fragments[0].seq_file.ino;
    assert!(fs.fragment_from_ino(ino).is_some());
}

#[test]
fn fasta_and_seq_inodes_are_distinct() {
    let f = fasta(">chr1\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_ne!(fs.fragments[0].fasta_file.ino, fs.fragments[0].seq_file.ino);
}

#[test]
fn all_fragment_inodes_are_unique() {
    let f = fasta(">chr1\nACGT\n>chr2\nGGGG\n>chr3\nTTTT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    let mut inos: Vec<u64> = fs
        .fragments
        .iter()
        .flat_map(|frag| [frag.fasta_file.ino, frag.seq_file.ino])
        .collect();
    let total = inos.len();
    inos.sort_unstable();
    inos.dedup();
    assert_eq!(inos.len(), total, "inode collision detected");
}

#[test]
fn label_no_name() {
    let f = fasta(">chr1\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(fs.fragments[0].label(), ">chr1\n");
}

#[test]
fn label_with_name() {
    let f = fasta(">chr1 Human chromosome 1\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(fs.fragments[0].label(), ">chr1 Human chromosome 1\n");
}

#[test]
fn fragment_name_parsed_correctly() {
    let f = fasta(">seq1 description text\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(fs.fragments[0].name.as_deref(), Some("description text"));
}

#[test]
fn data_size_ram_single_line() {
    // RAM cache stores pure sequence without newlines.
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(fs.fragments[0].data_size(), 8);
}

#[test]
fn data_size_ram_multiline() {
    // Newlines are stripped at read time; 8 bases regardless of line breaks.
    let f = fasta(">seq1\nACGT\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(fs.fragments[0].data_size(), 8);
}

#[test]
fn data_size_file_includes_newlines() {
    // File cache spans the raw byte region, newlines included.
    let f = fasta(">seq1\nACGT\nACGT\n");
    let fs = FustaFS::new(file_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(fs.fragments[0].data_size(), 10); // "ACGT\nACGT\n"
}

#[test]
fn chunk_ram_full_sequence() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].chunk(0, 8), b"ACGTACGT");
}

#[test]
fn chunk_ram_partial_with_offset() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].chunk(4, 4), b"ACGT");
}

#[test]
fn chunk_file_full_sequence() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(file_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].chunk(0, 8), b"ACGTACGT");
}

#[test]
fn chunk_file_partial_with_offset() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(file_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].chunk(2, 4), b"GTAC");
}

#[test]
fn chunk_file_second_fragment() {
    let f = fasta(">seq1\nAAAA\n>seq2\nGGGG\n");
    let fs = FustaFS::new(file_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragment_from_id("seq2").unwrap().chunk(0, 4), b"GGGG");
}

#[test]
fn chunk_mmap_full_sequence() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(mmap_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].chunk(0, 8), b"ACGTACGT");
}

#[test]
fn chunk_mmap_partial_with_offset() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(mmap_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].chunk(2, 4), b"GTAC");
}

#[test]
fn pure_chunk_ram_full_sequence() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].pure_chunk(0, 8), b"ACGTACGT");
}

#[test]
fn pure_chunk_ram_with_offset() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].pure_chunk(4, 4), b"ACGT");
}

#[test]
fn pure_chunk_mmap_strips_internal_newlines() {
    // Raw region is "ACGT\nACGT\n"; pure_chunk must return the 8 bases without \n.
    let f = fasta(">seq1\nACGT\nACGT\n");
    let fs = FustaFS::new(mmap_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].pure_chunk(0, 8), b"ACGTACGT");
}

#[test]
fn pure_chunk_mmap_offset_skips_past_newline() {
    let f = fasta(">seq1\nACGT\nACGT\n");
    let fs = FustaFS::new(mmap_settings(), f.path().to_str().unwrap()).unwrap();
    // Offset 4 in the newline-free view is the start of the second ACGT.
    assert_eq!(&*fs.fragments[0].pure_chunk(4, 4), b"ACGT");
}

#[test]
fn pure_chunk_file_first_fragment() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(file_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].pure_chunk(0, 8), b"ACGTACGT");
}

#[test]
fn pure_chunk_file_with_offset() {
    let f = fasta(">seq1\nACGTACGT\n");
    let fs = FustaFS::new(file_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].pure_chunk(4, 4), b"ACGT");
}

#[test]
fn pure_chunk_file_second_fragment() {
    let f = fasta(">seq1\nAAAA\n>seq2\nGGGG\n");
    let fs = FustaFS::new(file_settings(), f.path().to_str().unwrap()).unwrap();
    let frag2 = fs.fragment_from_id("seq2").unwrap();
    assert_eq!(&*frag2.pure_chunk(0, 4), b"GGGG");
}

#[test]
fn pure_chunk_file_strips_internal_newlines() {
    let f = fasta(">seq1\nACGT\nACGT\n");
    let fs = FustaFS::new(file_settings(), f.path().to_str().unwrap()).unwrap();
    assert_eq!(&*fs.fragments[0].pure_chunk(0, 8), b"ACGTACGT");
}

#[test]
fn fasta_ino_classified_correctly() {
    let f = fasta(">chr1\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    let ino = fs.fragments[0].fasta_file.ino;
    assert!(fs.is_fasta_file(ino));
    assert!(!fs.is_seq_file(ino));
    assert!(!fs.is_writeable(ino));
}

#[test]
fn seq_ino_classified_correctly() {
    let f = fasta(">chr1\nACGT\n");
    let fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    let ino = fs.fragments[0].seq_file.ino;
    assert!(fs.is_seq_file(ino));
    assert!(!fs.is_fasta_file(ino));
    assert!(fs.is_writeable(ino));
}

#[test]
fn labels_buffer_contains_all_headers() {
    let f = fasta(">chr1\nACGT\n>chr2 Human\nGGGG\n");
    let mut fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    fs.make_labels_buffer();
    let data = fs.get_file(LABELS_FILE).unwrap().data().to_vec();
    let text = String::from_utf8(data).unwrap();
    assert!(text.contains(">chr1\n"), "chr1 header missing");
    assert!(text.contains(">chr2 Human\n"), "chr2 header missing");
}

#[test]
fn labels_buffer_omits_sequence_data() {
    let f = fasta(">chr1\nACGT\n");
    let mut fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    fs.make_labels_buffer();
    let data = fs.get_file(LABELS_FILE).unwrap().data().to_vec();
    let text = String::from_utf8(data).unwrap();
    assert!(
        !text.contains("ACGT"),
        "sequence data leaked into labels buffer"
    );
}

#[test]
fn info_csv_uses_configured_separator() {
    let f = fasta(">chr1\nACGT\n");
    let mut fs = FustaFS::new(
        FustaSettings {
            csv_separator: "\t".to_string(),
            ..ram_settings()
        },
        f.path().to_str().unwrap(),
    )
    .unwrap();
    fs.make_info_csv_buffer();
    let data = fs.get_file(INFO_CSV_FILE).unwrap().data().to_vec();
    let text = String::from_utf8(data).unwrap();
    assert!(text.contains('\t'), "tab separator missing from CSV output");
}

#[test]
fn info_csv_lists_every_fragment() {
    let f = fasta(">chr1\nACGT\n>chr2\nGGGG\n");
    let mut fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    fs.make_info_csv_buffer();
    let data = fs.get_file(INFO_CSV_FILE).unwrap().data().to_vec();
    let text = String::from_utf8(data).unwrap();
    assert!(text.contains("chr1"));
    assert!(text.contains("chr2"));
}

#[test]
fn info_csv_has_header_row_with_id_and_length() {
    let f = fasta(">chr1\nACGT\n");
    let mut fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    fs.make_info_csv_buffer();
    let data = fs.get_file(INFO_CSV_FILE).unwrap().data().to_vec();
    let text = String::from_utf8(data).unwrap();
    let first = text.lines().next().unwrap();
    assert!(
        first.contains("id") && first.contains("length"),
        "unexpected CSV header: {:?}",
        first
    );
}

#[test]
fn info_csv_reports_correct_length() {
    let f = fasta(">seq1\nACGTACGT\n");
    let mut fs = FustaFS::new(ram_settings(), f.path().to_str().unwrap()).unwrap();
    fs.make_info_csv_buffer();
    let data = fs.get_file(INFO_CSV_FILE).unwrap().data().to_vec();
    let text = String::from_utf8(data).unwrap();
    // RAM cache stores pure sequence: 8 bases.
    assert!(text.contains("8"), "sequence length missing from CSV");
}
