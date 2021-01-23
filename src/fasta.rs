use std::io::prelude::*;
use std::io::{BufReader, Lines};

#[derive(Debug)]
pub struct Fragment {
    pub id: smartstring::SmartString<smartstring::LazyCompact>,
    pub name: Option<String>,
    pub pos: (usize, usize),
    pub len: usize,
}

pub struct FastaReader<T> {
    buffer_lines: Lines<BufReader<T>>,
    current_header: Option<String>,
    current_start: usize,
    current_offset: usize,
}

impl<T: Read> FastaReader<T> {
    pub fn new(file: T) -> FastaReader<T> {
        FastaReader {
            buffer_lines: BufReader::new(file).lines(),
            current_header: None,
            current_start: 0,
            current_offset: 0,
        }
    }
}

impl<T: Read> Iterator for FastaReader<T> {
    type Item = Fragment;

    fn next(&mut self) -> Option<Fragment> {
        while let Some(l) = self.buffer_lines.next() {
            let line = l.unwrap();
            let len = line.len() + 1;
            self.current_offset += len;

            if line.starts_with('>') {
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
                    };
                    self.current_header = Some(String::from(&line[1..]));

                    self.current_start = self.current_offset;

                    return Some(r);
                } else {
                    self.current_header = Some(String::from(&line[1..]));
                    self.current_start = self.current_offset;
                }
                continue;
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
            };
            self.current_header = None;
            return Some(r);
        }

        None
    }
}
