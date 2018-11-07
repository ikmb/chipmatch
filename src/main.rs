extern crate zip;

#[macro_use]
extern crate clap;

use std::cmp::{Ordering, PartialEq, PartialOrd};
use std::collections::BinaryHeap;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Result, Write};

use clap::{App, Arg};

/// When are f32's considered equal?
/// Needed to establish some flavor of total ordering on floats
const F32_EPSILON: f32 = 0.00001;
const EXTRACT_BUFFER_SIZE: usize = 1024 * 1024;

/// The result structure for a single strand file match
struct MatchResult {
    /// Name of the actual strand file within the ZIP archive
    pub name: String,

    /// Ratio of how many BIM file entries have successfully matched against the strand file
    pub name_match_rate: f32,
    pub name_pos_match_rate: f32,
    pub strand_match_rate: f32,
    pub plus_match_rate: f32,
    pub atcg_match_rate: f32,
}

#[derive(Clone)]
struct SourceEntry {
    pub name: String,
    pub chromosome: u64,
    pub position: u64,
    pub alleles: (char, char),
}

#[derive(Clone)]
struct VariantEntry {
    pub chromosome: u64,
    pub position: u64,
    pub alleles: (char, char),
    pub strand: char,
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
        let my_score = self.name_match_rate * 4000.0
            + self.name_pos_match_rate * 3000.0
            + self.strand_match_rate * 2000.0
            + self.plus_match_rate;

        let other_score = other.name_match_rate * 4000.0
            + other.name_pos_match_rate * 3000.0
            + other.strand_match_rate * 2000.0
            + other.plus_match_rate;

        if my_score < other_score {
            Some(Ordering::Less)
        } else if my_score > other_score {
            Some(Ordering::Greater)
        } else if (my_score - other_score).abs() < F32_EPSILON {
            Some(Ordering::Equal)
        } else {
            eprintln!("score: {} {}", my_score, other_score);
            panic!("Aaaah!");
        }
    }
}

impl PartialEq for MatchResult {
    fn eq(&self, other: &MatchResult) -> bool {
        self.cmp(&other) == Ordering::Equal
    }
}

fn chromosome_to_number(s: &str) -> u64 {
    let n = if let Ok(num) = u64::from_str_radix(&s, 10) {
        num
    } else {
        match s {
            "X" => 23,
            "Y" => 24,
            "XY" => 25,
            "M" | "MT" => 26,
            _ => 0,
        }
    };

    n
}

// Reads a list of variants from the BIM file
fn read_bim(filename: &str) -> Result<Vec<SourceEntry>> {
    let file = File::open(filename)?;
    let mut names: Vec<SourceEntry> = Vec::new();
    let mut entry = SourceEntry {
        name: String::from(""),
        chromosome: 0,
        position: 0,
        alleles: ('X', 'X'),
    };

    // Data syntax:
    // <some value> <variant-name*> <some other value> <position*> <some other value> ...
    //              ^^^^^^^^^^^^^^^ String             ^^^^^^^^^^^^ u64
    for line in BufReader::new(file).lines() {
        let line = line?;
        let mut l = line.split_whitespace();
        entry.chromosome = chromosome_to_number(l.next().unwrap());
        entry.name = l.next().unwrap().to_string();
        entry.position = u64::from_str_radix(l.nth(1).unwrap(), 10).unwrap();
        entry.alleles.0 = l.next().unwrap().chars().next().unwrap();
        entry.alleles.1 = l.next().unwrap().chars().next().unwrap();

        names.push(entry.clone());
    }
    Ok(names)
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
fn read_variants_from_zip(filename: &str) -> Result<(String, HashMap<String, VariantEntry>)> {
    let mut zip = zip::ZipArchive::new(File::open(filename)?)?;
    let mut variants = HashMap::new();
    let mut strand_file_name: String = String::new();

    let mut var: VariantEntry = VariantEntry {
        chromosome: 0,
        position: 0,
        alleles: ('X', 'X'),
        strand: 'X',
    };

    for i in 0..zip.len() {
        let mut file = zip.by_index(i).unwrap();
        if file.name().ends_with(".strand") {
            strand_file_name = file.name().to_string();
            for line in BufReader::new(file).lines() {
                let line = line?;
                let mut l = line.split_whitespace();

                let name = l.next().unwrap().to_string();
                var.chromosome = chromosome_to_number(l.next().unwrap());
                var.position = u64::from_str_radix(l.next().unwrap(), 10).unwrap_or(0);
                l.next().unwrap_or("*");
                var.strand = l.next().unwrap().chars().next().unwrap();
                //                println!("{} {} {}", var.chromosome, var.position, var.strand);

                // Not all strand files carry allele information
                let mut alleles = l.next().unwrap_or("XX").chars();
                var.alleles.0 = alleles.next().unwrap();
                var.alleles.1 = alleles.next().unwrap();
                variants.insert(name, var.clone());
            }
            // We don't need more than one strand file
            break;
        }
    }

    Ok((strand_file_name, variants))
}

enum AlleleMatch {
    Original,
    Plus,
    ATCG,
    Mismatch,
}

fn match_set(left: (char, char), right: (char, char)) -> bool {
    (left.0 == right.0 || left.0 == right.1) && (left.1 == right.0 || left.1 == right.1)
}

fn flip_alleles(i: (char, char)) -> (char, char) {
    let mut res: (char, char) = (' ', ' ');

    res.0 = match i.0 {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => 'X',
    };
    res.1 = match i.1 {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => 'X',
    };
    res
}

fn match_alleles(left: (char, char), right: (char, char), strand: char) -> AlleleMatch {
    if match_set(left, flip_alleles(left)) && match_set(right, flip_alleles(right)) {
        AlleleMatch::ATCG
    } else if match_set(left, right) {
        AlleleMatch::Original
    } else if match_set(left, flip_alleles(right)) && strand == '-' {
        AlleleMatch::Plus
    } else {
        AlleleMatch::Mismatch
    }
}

fn mzinf(score: f32) -> f32 {
    if score.is_finite() {
        score
    } else {
        0.0
    }
}

// Match a BIM file against a strand file
fn match_bim(
    bim: &[SourceEntry],
    name: &str,
    variants: &HashMap<String, VariantEntry>,
) -> MatchResult {
    let mut res = MatchResult {
        name: name.to_string(),
        name_match_rate: 0.0,
        name_pos_match_rate: 0.0,
        strand_match_rate: 0.0,
        plus_match_rate: 0.0,
        atcg_match_rate: 0.0,
    };

    let mut name_matches = 0;
    let mut name_pos_matches = 0;

    let mut strand_matches = 0;
    let mut plus_matches = 0;
    let mut atcg_matches = 0;

    for bimentry in bim {
        if let Some(strand) = variants.get(&bimentry.name) {
            name_matches += 1;

            if (strand.position == bimentry.position) && (strand.chromosome == bimentry.chromosome)
            {
                name_pos_matches += 1;

                match match_alleles(bimentry.alleles, strand.alleles, strand.strand) {
                    AlleleMatch::ATCG => atcg_matches += 1,
                    AlleleMatch::Original => strand_matches += 1,
                    AlleleMatch::Plus => plus_matches += 1,
                    AlleleMatch::Mismatch => {}
                }
            }
        }
    }

    res.name = name.to_string();
    res.name_match_rate = name_matches as f32 / bim.len() as f32;
    res.name_pos_match_rate = name_pos_matches as f32 / name_matches as f32;
    res.strand_match_rate = strand_matches as f32 / (name_pos_matches - atcg_matches) as f32;
    res.plus_match_rate = plus_matches as f32 / (name_pos_matches - atcg_matches) as f32;
    res.atcg_match_rate = atcg_matches as f32 / name_pos_matches as f32;

    res.name_match_rate = mzinf(res.name_match_rate);
    res.name_pos_match_rate = mzinf(res.name_pos_match_rate);
    res.strand_match_rate = mzinf(res.strand_match_rate);
    res.plus_match_rate = mzinf(res.plus_match_rate);
    res.atcg_match_rate = mzinf(res.atcg_match_rate);

    res
}

// Extracts a strand file from the given ZIP archive
// and dumps it to the current working directory by
// its original name
fn extract_strand(zipfile: &str) -> Result<()> {
    let mut zip = zip::ZipArchive::new(File::open(zipfile)?)?;

    let mut buffer = vec![0; EXTRACT_BUFFER_SIZE];

    for i in 0..zip.len() {
        let mut file = zip.by_index(i).unwrap();
        if file.name().ends_with(".strand") {
            let mut target = File::create(file.name())?;
            while let Ok(size) = file.read(&mut buffer) {
                if size == 0 {
                    break;
                }
                target.write_all(&buffer[0..size])?;
            }
            break;
        }
    }

    Ok(())
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
                .help("Directory containing Will Rayner's strand archives"),
        )
        .arg(
            Arg::with_name("verbose")
                .short("v")
                .long("verbose")
                .help("Be verbose and print progress"),
        )
        .arg(
            Arg::with_name("extract")
                .short("e")
                .long("extract")
                .value_name("N")
                .help("Extract the N most promising strand files to the local working directory")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .takes_value(true)
                .value_name("FILE")
                .help("Write result table to FILE instead of stdout"),
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
    let mut strandmap: HashMap<String, String> = HashMap::new();

    for z in ziplist {
        if verbose > 0 {
            println!("Scanning {}", z);
        }

        let (name, vars) = read_variants_from_zip(&z)?;
        strandmap.insert(name.to_string(), z);
        let res = match_bim(&bim, &name, &vars);
        results.push(res);
    }

    let mut extract_strands =
        u64::from_str_radix(matches.value_of("extract").unwrap_or("0"), 10).unwrap();

    // Set output target
    let mut out_writer: Box<Write> = match matches.value_of("output") {
        Some(filename) => Box::new(File::create(&filename)?),
        None => Box::new(io::stdout()),
    };

    writeln!(&mut out_writer, "strand\tname_match_rate\tpos_match_rate\toriginal_match_rate\tplus_match_rate\tatcg_match_rate")?;

    while let Some(res) = results.pop() {
        if extract_strands > 0 {
            if verbose > 0 {
                println!("Extracting {} from {}", &res.name, &strandmap[&res.name]);
            }
            extract_strand(&strandmap[&res.name])?;
            extract_strands -= 1;
        }

        writeln!(
            &mut out_writer,
            "{}\t{}\t{}\t{}\t{}\t{}",
            res.name,
            res.name_match_rate,
            res.name_pos_match_rate,
            res.strand_match_rate,
            res.strand_match_rate + res.plus_match_rate,
            res.atcg_match_rate
        )?;
    }

    Ok(())
}
