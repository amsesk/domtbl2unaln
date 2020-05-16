extern crate bio;
use bio::utils::Text;
use clap::{App, Arg};
use csv::WriterBuilder;
use lazy_static::lazy_static;
use multimap::MultiMap;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::fmt;
use std::fs::File;

use std::io::{BufRead, Write};
use std::io::{BufReader, BufWriter};
use std::path::Path;

pub fn main() -> Result<(), std::io::Error> {
    let args = App::new("domtbl-reader")
        .version("0.1")
        .author("Kevin Amses")
        .about("Reads and parses domtbl from hmmsearch.")
        .arg(
            Arg::with_name("domtbls")
                .long("domtbls")
                .value_name("DIR")
                .about("Path to directory containing domtbls.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("proteins")
                .long("proteins")
                .value_name("DIR")
                .about("A protein fasta including predictions for ALL taxa, with samtools faidx.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("outdir")
                .long("outdir")
                .value_name("DIR")
                .about("Path to existing directory to dump unaligned fastas.")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let mut hitmap = HashMap::new();

    // Unwrap on None here is impossible because
    // clap will have already killed the process
    // since a required argument is missing
    let domtbls = Path::new(args.value_of("domtbls").unwrap());
    let proteins = Path::new(args.value_of("proteins").unwrap());
    let outdir = Path::new(args.value_of("outdir").unwrap())
        .canonicalize()
        .unwrap();

    let paths = file_list(&domtbls, "domtbl").unwrap();
    let hits: Vec<(&std::path::PathBuf, Hits)> = paths
        .iter()
        .map(|p| (p, parse_and_filter(p).unwrap()))
        .collect();

    for h in hits {
        hitmap.insert(h.0, h.1);
    }

    // Write marker recovery rates by taxon
    let mut wtr = WriterBuilder::new()
        .delimiter(b'\t')
        .from_path("/home/aimzez/DATA/phylogeny/marker_recovery.csv")?;
    for (k, v) in hitmap.iter() {
        let fname_repr = k.to_str().to_owned().unwrap();
        wtr.write_record(&[
            fname_repr,
            &format!("{}", v.inner.len()),
            &format!("{:.2}", v.perc_duplicated()),
        ])?;
    }
    wtr.flush()?;

    let odb10_marker_cutoffs = parse_cutoffs(ODB10_CUTOFFS);

    let mut index_reader = bio::io::fasta::IndexedReader::from_file(&proteins).unwrap();
    for m in odb10_marker_cutoffs.keys() {
        if !has_enough_occupants(&m, &hitmap, 0.75, 140.0) {
            continue;
        }
        let mut marker_fname = String::from(*m);
        marker_fname.push_str(".fasta");
        let outpath = outdir.join(Path::new(&marker_fname));
        let unaln = File::create(&outpath)?;
        let mut unaln_writer = BufWriter::new(&unaln);
        for (_k, h) in &hitmap {
            match h.hits_for(&m) {
                Some(hits) => {
                    for hit in hits {
                        let mut text = Text::new();
                        match index_reader.fetch_all(&hit.target) {
                            Ok(()) => {
                                index_reader.read(&mut text).unwrap();
                                write!(
                                    &mut unaln_writer,
                                    ">{}\n{}\n",
                                    &hit.target,
                                    String::from_utf8(text).unwrap()
                                )
                                .unwrap();
                            }
                            Err(_e) => {
                                println!("Did not find contig: {}", &hit.target);
                            }
                        }
                    }
                }
                None => (),
            }
        }
    }

    Ok(())
}

static ODB10_CUTOFFS: &str = include_str!("../lib/odb10_scores_cutoff");

lazy_static! {
    static ref DOMTBL_COLUMNS: HashMap<&'static str, usize> = {
        let mut map = HashMap::new();
        map.insert("target_name", 0);
        map.insert("target_accession", 1);
        map.insert("target_tlen", 2);
        map.insert("query_name", 3);
        map.insert("query_accession", 4);
        map.insert("qlen", 5);
        map.insert("fs_evalue", 6);
        map.insert("fs_score", 7);
        map.insert("fs_bias", 8);
        map.insert("dom_id", 9);
        map.insert("dom_total", 10);
        map.insert("dom_cevalue", 11);
        map.insert("dom_ievalue", 12);
        map.insert("dom_score", 13);
        map.insert("dom_bias", 14);
        map.insert("dom_id", 15);
        map.insert("hmm_from", 16);
        map.insert("hmm_to", 17);
        map.insert("ali_from", 18);
        map.insert("ali_to", 19);
        map.insert("env_from", 20);
        map.insert("env_to", 21);
        map.insert("acc", 22);
        map.insert("target_description", 23);

        map
    };
}

// Struct to hold a vector of
// hits for each key
// the keys are the markers here
#[derive(Clone)]
pub struct Hits {
    inner: MultiMap<String, Hit>,
    n_markers: usize,
}
impl Hits {
    pub fn from_file<P: AsRef<Path>>(_path: P) {
        ()
    }

    pub fn hits_for(&self, key: &str) -> Option<&Vec<Hit>> {
        self.inner.get_vec(key)
    }

    pub fn perc_duplicated(&self) -> f64 {
        let mut count: f64 = 0.0;
        for (_k, m) in self.inner.iter_all() {
            match m.len() {
                0 => panic!("This sholdn't have happened. Zero-length vector in Hits.inner"),
                1 => (),
                _ => {
                    count = count + 1.0;
                }
            }
        }
        let perc_duplicated = count / (self.inner.len() as f64);
        perc_duplicated
    }
}

impl TryFrom<MultiMap<String, Hit>> for Hits {
    type Error = &'static str;

    fn try_from(mmap: MultiMap<String, Hit>) -> Result<Hits, Self::Error> {
        Ok(Hits {
            inner: mmap.to_owned(),
            n_markers: mmap.len(),
        })
    }
}

// Struct for each hit that holds
// important columns as well as a
// vector of all of the columns
#[derive(Clone)]
pub struct Hit {
    pub target: String,
    pub query: String,
    pub fs_evalue: f64,
    pub fs_score: f64,
    pub inner: Vec<String>,
}

impl Hit {
    pub fn new(spl: Vec<String>) -> Hit {
        Hit {
            target: spl[DOMTBL_COLUMNS["target_name"]].to_owned(),
            query: spl[DOMTBL_COLUMNS["query_name"]].to_owned(),
            fs_evalue: spl[DOMTBL_COLUMNS["fs_evalue"]]
                .to_owned()
                .parse::<f64>()
                .expect("Could not convert fs_evalue to f64"),
            fs_score: spl[DOMTBL_COLUMNS["fs_score"]]
                .to_owned()
                .parse::<f64>()
                .expect("Could not convert fs_score to f64"),
            inner: spl,
        }
    }
}

impl fmt::Debug for Hit {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "{}", self.target)
    }
}

pub fn file_list<'a>(
    path: &'a Path,
    ext: &'static str,
) -> Result<Vec<std::path::PathBuf>, std::io::Error> {
    let ls = path.read_dir()?;
    let lsfilt = ls
        .map(|p| p.unwrap().path())
        .filter(|p| !p.is_dir())
        .filter(|p| p.extension().unwrap() == ext)
        .map(|p| p.canonicalize().unwrap())
        .collect();
    Ok(lsfilt)
}

fn has_enough_occupants(
    marker: &str,
    hitmap: &HashMap<&std::path::PathBuf, Hits>,
    proportion: f64,
    n_taxa: f64,
) -> bool {
    let mut i = 0.0;
    for (_k, h) in hitmap {
        match h.hits_for(&marker) {
            Some(_hits) => i += 1.0,
            None => (),
        }
    }
    let mark_occ = i / n_taxa;
    println!("{}", mark_occ);
    if mark_occ <= proportion {
        return false;
    }
    true
}

pub fn parse_cutoffs(cutoffs: &'static str) -> HashMap<&'static str, f64> {
    let mut cutoffs_map = HashMap::new();
    for line in cutoffs.split('\n').collect::<Vec<&'static str>>() {
        if line.is_empty() {
            continue;
        }
        let spl = line.split('\t').collect::<Vec<&'static str>>();
        cutoffs_map.insert(spl[0], spl[1].parse::<f64>().unwrap());
    }
    cutoffs_map
}

//pub fn parse_and_filter() -> Result<MultiMap<String, Hit>, std::io::Error> {
pub fn parse_and_filter<P>(domtbl_path: P) -> Result<Hits, std::io::Error>
where
    P: AsRef<Path>,
{
    let domtbl = File::open(domtbl_path)?;
    let domtbl = BufReader::new(domtbl);

    // Iterate over lines in domtbl and add values to MultiMap with query_name as key
    let mut markerhits = MultiMap::new();
    for line in domtbl
        .lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            l.split_ascii_whitespace()
                .map(|p| p.to_owned())
                .collect::<Vec<String>>()
        })
        .into_iter()
    {
        let thishit = Hit::new(line);
        markerhits.insert(String::from(&thishit.query), thishit);
    }

    let odb10_cutoffs = parse_cutoffs(ODB10_CUTOFFS);

    let mut filt = filter_by_score(&odb10_cutoffs, markerhits);

    dedup_hits(&mut filt);

    Ok(Hits::try_from(filt).unwrap())
}

pub fn filter_by_score(
    cutoffs: &HashMap<&'static str, f64>,
    hits: MultiMap<String, Hit>,
) -> MultiMap<String, Hit> {
    let mut new = MultiMap::new();
    for (key, values) in hits.iter_all() {
        for val in values.iter() {
            if val.fs_score >= *cutoffs.get(&key as &str).unwrap() {
                new.insert(String::from(key), val.to_owned());
            }
        }
    }
    new
}

pub fn dedup_hits(hits: &mut MultiMap<String, Hit>) {
    for (_k, v) in hits.iter_all_mut() {
        v.sort_by(|a, b| a.target.cmp(&b.target));
        v.dedup_by(|a, b| {
            (a.target == b.target)
                && (a.query == b.query)
                && (a.fs_evalue == b.fs_evalue)
                && (a.fs_score == b.fs_score)
        });
    }
}
