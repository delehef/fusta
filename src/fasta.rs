use std::io::prelude::*;
use std::io::{BufReader, Lines};

#[derive(Debug)]
pub struct Fragment {
    pub id: smartstring::SmartString<smartstring::LazyCompact>,
    pub name: Option<String>,
    pub pos: (usize, usize),
    pub len: usize,
    pub seq: Option<Vec<u8>>,
}

pub struct FastaReader<T> {
    buffer_lines: Lines<BufReader<T>>,
    current_header: Option<String>,
    current_start: usize,
    current_offset: usize,

    with_seq: bool,
}

impl<T: Read> FastaReader<T> {
    pub fn new(file: T, with_seq: bool) -> FastaReader<T> {
        FastaReader {
            buffer_lines: BufReader::new(file).lines(),
            current_header: None,
            current_start: 0,
            current_offset: 0,

            with_seq,
        }
    }
}

impl<T: Read> Iterator for FastaReader<T> {
    type Item = Fragment;

    fn next(&mut self) -> Option<Fragment> {
        let mut current_seq: Vec<u8> = Vec::new();
        for l in self.buffer_lines.by_ref() {
            let line = l.unwrap();
            let len = line.len() + 1;
            self.current_offset += len;

            if let Some(name) = line.strip_prefix('>') {
                if let Some(ref current_header) = self.current_header {
                    let split = current_header.split(' ').collect::<Vec<_>>();
                    let r = Fragment {
                        id: split[0].into(),
                        name: if split.len() > 1 {
                            Some(split[1..].join(" "))
                        } else {
                            None
                        },

                        pos: (self.current_start, self.current_offset - len),
                        len: self.current_offset - self.current_start - len,
                        seq: if self.with_seq {
                            Some(current_seq)
                        } else {
                            None
                        },
                    };
                    self.current_header = Some(String::from(name));

                    self.current_start = self.current_offset;

                    return Some(r);
                } else {
                    self.current_header = Some(String::from(name));
                    self.current_start = self.current_offset;
                }
                continue;
            } else if self.with_seq {
                current_seq.extend(line.trim_end().as_bytes());
            }
        }

        if let Some(ref current_header) = self.current_header {
            let split = current_header.split(' ').collect::<Vec<_>>();
            let r = Fragment {
                id: split[0].into(),
                name: if split.len() > 1 {
                    Some(split[1..].join(" "))
                } else {
                    None
                },
                pos: (self.current_start, self.current_offset),
                len: self.current_offset - self.current_start,
                seq: if self.with_seq {
                    Some(current_seq)
                } else {
                    None
                },
            };
            self.current_header = None;
            return Some(r);
        }

        None
    }
}
