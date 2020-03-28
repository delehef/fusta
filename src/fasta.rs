use std::io::{BufReader, Lines};
use std::io::prelude::*;
use std::fs::File;

pub struct Fragment {
    pub id: String,
    pub name: Option<String>,
    pub sequence: Vec<u8>
}

pub struct FastaReader<T> {
    buffer_lines: Lines<BufReader<T>>,
    current_header: Option<String>,
    current_sequence: String
}

impl<T: Read> FastaReader<T> {
    fn new(file: T) -> FastaReader<T> {
        FastaReader {
            buffer_lines: BufReader::new(file).lines(),
            current_header: None,
            current_sequence: String::new()
        }
    }
}

impl<T: Read> Iterator for FastaReader<T> {
    type Item = Fragment;

    fn next(&mut self) -> Option<Fragment> {
        while let Some(l) = self.buffer_lines.next() {
            let line = l.unwrap();
            if line.starts_with(">") {
                if let Some(ref current_header) = self.current_header {
                    let split = current_header.split(" ").collect::<Vec<_>>();
                    let r = Fragment {
                        id: split[0].to_string(),
                        name: if split.len() > 1 { Some(split[1..].join(" ")) } else { None },
                        sequence: self.current_sequence.as_bytes().to_vec(),
                    };
                    self.current_header = Some(String::from(&line[1..]));
                    self.current_sequence.clear();
                    return Some(r);
                } else {
                    self.current_header = Some(String::from(&line[1..]));
                    self.current_sequence.clear();
                }
                continue;
            }
            self.current_sequence.push_str(line.trim());
        }

        if let Some(ref current_header) = self.current_header {
            let split = current_header.split(" ").collect::<Vec<_>>();
            let r = Fragment {
                id: split[0].to_string(),
                name: if split.len() > 1 { Some(split[1..].join(" ")) } else { None },
                sequence: self.current_sequence.as_bytes().to_vec(),
            };
            self.current_header = None;
            self.current_sequence.clear();
            self.current_sequence.shrink_to_fit();
            return Some(r);
        }

        None
    }
}


pub fn from_file(filename: &str) -> Result<Vec<Fragment>, &str> {
    let file = File::open(filename).or_else(|_| Err("Unable to open FASTA file")).unwrap();
    Ok(FastaReader::new(file).collect::<Vec<_>>())
}
