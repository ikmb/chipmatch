extern crate zip;

#[macro_use]
extern crate clap;

use std::cmp::{Ordering, PartialEq, PartialOrd};
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Result};

use clap::{App, Arg};

/// When are f32's considered equal?
/// Needed to establish some flavor of total ordering on floats
const F32_EPSILON: f32 = 0.00001;

/// The result structure for a single strand file match
struct MatchResult {
    /// Name of the actual strand file within the ZIP archive
    pub name: String,

    /// Ratio of how many BIM file entries have successfully matched against the strand file
    pub match_rate: f32,

    /// Ratio of BIM file matches vs. total number of variants in strand file
    pub completeness: f32,
}

// Define a total order with the help of a partial order on the MatchResult
// so we can sort the MatchResults later by match_rate and completeness
impl Eq for MatchResult {}

impl Ord for MatchResult {
    fn cmp(&self, other: &MatchResult) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

// Define a partial order on the MatchResult
impl PartialOrd for MatchResult {
    fn partial_cmp(&self, other: &MatchResult) -> Option<Ordering> {
        if self.match_rate < other.match_rate {
            return Some(Ordering::Less);
        } else if self.match_rate > other.match_rate {
            return Some(Ordering::Greater);
        } else if (self.match_rate - other.match_rate).abs() < F32_EPSILON {
            if self.completeness < other.completeness {
                return Some(Ordering::Less);
            } else if self.completeness > other.completeness {
                return Some(Ordering::Greater);
            } else if (self.completeness - other.completeness).abs() < F32_EPSILON {
                return Some(Ordering::Equal);
            } else {
                panic!("Completeness is NaN or inf");
            }
        } else {
            panic!("Match rate is NaN or inf");
        }
    }
}

impl PartialEq for MatchResult {
    fn eq(&self, other: &MatchResult) -> bool {
        (self.match_rate - other.match_rate).abs() < F32_EPSILON
            && (self.completeness - other.completeness).abs() < F32_EPSILON
    }
}

// Reads a list of variants from the BIM file
fn read_bim(filename: &str) -> Result<Vec<(String, u64)>> {
    let file = File::open(filename)?;
    let mut names: Vec<(String, u64)> = Vec::new();

    // Data syntax:
    // <some value> <variant-name*> <some other value> <position*> <some other value> ...
    //              ^^^^^^^^^^^^^^^ String             ^^^^^^^^^^^^ u64
    for line in BufReader::new(file).lines() {
        let line = line?;
        let mut l = line.split_whitespace();
        let v = l.nth(1).unwrap();
        let pos = u64::from_str_radix(l.nth(1).unwrap(), 10).unwrap();
        names.push((v.to_string(), pos));
    }
    return Ok(names);
}

// Read a list of ZIP files from the given directory.
// Each ZIP file name has the directory name prepended
fn get_zip_list(dirname: &str) -> Result<Vec<String>> {
    let entries = fs::read_dir(dirname)?;
    let mut names: Vec<String> = Vec::new();
    for entry in entries {
        let f = entry?;
        if let Ok(name) = f.file_name().into_string() {
            if name.ends_with(".zip") {
                names.push(format!("{}/{}", dirname, name).to_string());
            }
        }
    }
    Ok(names)
}

/// Find the strand file within a ZIP archive and extract the variant/position pairs
fn read_variants_from_zip(filename: &str) -> Result<(String, HashMap<String, u64>)> {
    let mut zip = zip::ZipArchive::new(File::open(filename)?)?;
    let mut variants: HashMap<String, u64> = HashMap::new();
    let mut strand_file_name: String = String::new();

    for i in 0..zip.len() {
        let mut file = zip.by_index(i).unwrap();
        if file.name().ends_with(".strand") {
            strand_file_name = file.name().to_string();
            for line in BufReader::new(file).lines() {
                let line = line?;
                let mut l = line.split_whitespace();
                let v = l.next().unwrap().to_string();
                let pos = u64::from_str_radix(l.nth(1).unwrap(), 10).unwrap_or(0);
                variants.insert(v, pos);
            }
            // We don't need more than one strand file
            break;
        }
    }

    Ok((strand_file_name, variants))
}

fn match_bim(bim: &[(String, u64)], name: &str, variants: HashMap<String, u64>) -> MatchResult {
    let mut res = MatchResult {
        name: name.to_string(),
        match_rate: 0.0,
        completeness: 0.0,
    };

    let mut match_count = 0;

    for variant in bim {
        if let Some(pos) = variants.get(&variant.0 as &str) {
            if *pos == variant.1 {
                match_count += 1;
            }
        }
    }

    res.name = name.to_string();
    res.match_rate = match_count as f32 / bim.len() as f32;
    res.completeness = match_count as f32 / variants.len() as f32;

    if !res.match_rate.is_finite() {
        res.match_rate = 0.0;
    }

    if !res.completeness.is_finite() {
        res.completeness = 0.0;
    }

    res
}

fn main() -> Result<()> {
    let matches = App::new(crate_name!())
        .version(crate_version!())
        .author(crate_authors!())
        .about(crate_description!())
        .arg(
            Arg::with_name("bim")
                .value_name("FILE")
                .takes_value(true)
                .required(true)
                .help("PLINK .bim file to guess the chip type for"),
        )
        .arg(
            Arg::with_name("strandfolder")
                .value_name("DIR")
                .takes_value(true)
                .required(true)
                .help("Directory containing Will Rainer's strand archives"),
        )
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .help("Be verbose and print progress"),
        )
        .get_matches();

    let verbose = matches.occurrences_of("verbose");

    if verbose > 0 {
        println!("Reading BIM file...");
    }
    let bim = read_bim(matches.value_of("bim").unwrap())?;
    if verbose > 0 {
        println!("{} variants loaded.", bim.len());
    }
    let ziplist = get_zip_list(matches.value_of("strandfolder").unwrap())?;

    let mut results = BinaryHeap::new();

    for z in ziplist {
        if verbose > 0 {
            println!("Scanning {}", z);
        }

        let (name, vars) = read_variants_from_zip(&z)?;
        let res = match_bim(&bim, &name, vars);
        results.push(res);
    }

    while let Some(res) = results.pop() {
        println!("{}\t{}\t{}", res.name, res.match_rate, res.completeness);
    }

    Ok(())
}
