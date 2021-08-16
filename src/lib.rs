pub mod hits;
use crate::hits::{Hit, Hits};
use csv::WriterBuilder;
use lazy_static::lazy_static;
use multimap::MultiMap;
use std::collections::HashMap;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::BufRead;
use std::io::BufReader;
use std::path::{Path, PathBuf};

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

pub fn has_enough_occupants(
    marker: &str,
    all_hits: &Vec<Hits>,
    proportion: f64,
    n_taxa: f64,
) -> bool {
    let mut i = 0.0;
    for hits in all_hits {
        match hits.hits_for(&marker) {
            Some(_hits) => i += 1.0,
            None => (),
        }
    }
    let mark_occ = i / n_taxa;
    //println!("{}", mark_occ);
    if mark_occ < proportion {
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
pub fn parse_and_filter(
    domtbl_path: PathBuf,
    busco_filter: bool,
    cutoffs: &'static str,
    outdir: &Path,
) -> Result<Hits, std::io::Error> {
    println!("Parsing {}...", &domtbl_path.to_str().unwrap());
    let domtbl = File::open(&domtbl_path)?;
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

    let odb10_cutoffs = parse_cutoffs(cutoffs);

    if busco_filter {
        markerhits = filter_by_score(&odb10_cutoffs, markerhits);
    }

    calculate_aln_length(&markerhits, outdir);

    dedup_hits(&mut markerhits);

    let ret = Hits::new(
        domtbl_path.into_os_string().into_string().unwrap(),
        markerhits,
    );
    Ok(ret)
}

pub fn calculate_aln_length(markerhits: &MultiMap<String, Hit>, outdir: &Path) {
    let mut aln_lengths = MultiMap::new();
    for (marker, hits) in markerhits.iter_all() {
        for hit in hits {
            aln_lengths.insert((marker.to_owned(), hit.target.to_owned()), hit.aln_len)
        }
    }
    for (_key, lens) in aln_lengths.iter_all_mut() {
        let sum: i64 = lens.iter().sum();
        *lens = vec![sum];
    }
    println!("{:?}", outdir);
    let handle = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(false)
        .append(true)
        .open(outdir.join(Path::new("aln_length.csv")))
        .unwrap();

    let mut wtr = WriterBuilder::new().delimiter(b'\t').from_writer(handle);
    for (key, len) in aln_lengths.iter() {
        wtr.write_record(&[
            &format!("{}", key.0),
            &format!("{}", key.1),
            &format!("{}", len),
        ])
        .unwrap();
    }
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

pub static FUNGI_ODB10_CUTOFFS: &str = include_str!("../lib/fungi_odb10_scores_cutoff");
pub static MOLLICUTES_ODB10_CUTOFFS: &str = include_str!("../lib/mollicutes_odb10_scores_cutoff");
pub static BURKHOLDERIALES_ODB10_CUTOFFS: &str =
    include_str!("../lib/burkholderiales_odb10_scores_cutoff");
pub static JGI434_ODB10_CUTOFFS: &str = include_str!("../lib/jgi434_scores_cutoff");
pub static MUCORALES_ODB10_CUTOFFS: &str = include_str!("../lib/mucorales_odb10_scores_cutoff");

lazy_static! {
    static ref DOMTBL_COLUMNS: HashMap<&'static str, usize> = {
        let mut map = HashMap::new();
        map.insert("target_name", 0);
        map.insert("target_accession", 1);
        map.insert("tlen", 2);
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
        map.insert("hmm_from", 15);
        map.insert("hmm_to", 16);
        map.insert("ali_from", 17);
        map.insert("ali_to", 18);
        map.insert("env_from", 19);
        map.insert("env_to", 20);
        map.insert("acc", 21);
        map.insert("target_description", 22);

        map
    };
}
