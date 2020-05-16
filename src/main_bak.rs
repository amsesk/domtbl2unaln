extern crate bio;
use clap::{App, Arg};
use lazy_static::lazy_static;
use multimap::MultiMap;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

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

fn main() -> Result<(), std::io::Error> {
    let _domtbl_columns = [
        "target_name",
        "target_accession",
        "tlen",
        "query_name",
        "query_accession",
        "qlen",
    ];

    // Take command-line argument for --domtbl
    let matches = App::new("domtbl-reader")
        .version("0.1")
        .author("Kevin Amses")
        .about("Reads and parses domtbl from hmmsearch.")
        .arg(
            Arg::with_name("domtbl")
                .long("domtbl")
                .value_name("FILE")
                .about("Path to domtbl.")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let domtbl = File::open(matches.value_of("domtbl").unwrap())?;
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
        markerhits.insert(
            String::from(&line[DOMTBL_COLUMNS["query_name"]]),
            (
                //Only adding target_name, fs_evalue, and fs_score here
                String::from(&line[DOMTBL_COLUMNS["target_name"]]),
                String::from(&line[DOMTBL_COLUMNS["fs_evalue"]])
                    .parse::<f64>()
                    .unwrap(),
                String::from(&line[DOMTBL_COLUMNS["fs_score"]])
                    .parse::<f64>()
                    .unwrap(),
            ),
        )
    }

    let odb10_cutoffs = parse_cutoffs(ODB10_CUTOFFS);

    let mut filt = filter_by_score(&odb10_cutoffs, markerhits);
    for (_k, v) in filt.iter_all_mut() {
        v.sort_by(|(a, _b, _c), (d, _e, _f)| a.partial_cmp(d).unwrap());
        v.dedup();
    }
    println!("#marker\tprotein\tfs_evalue\tfs_score\tdup");
    for (k, v) in filt.iter_all() {
        for t in v {
            println!("{}\t{}\t{:e}\t{}\t{}", k, t.0, t.1, t.2, v.len());
        }
    }

    Ok(())
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

pub fn filter_by_score(
    cutoffs: &HashMap<&'static str, f64>,
    hits: MultiMap<String, (String, f64, f64)>,
) -> MultiMap<String, (String, f64, f64)> {
    let mut new = MultiMap::new();
    for (key, values) in hits.iter_all() {
        for val in values {
            if val.2 >= *cutoffs.get(&key as &str).unwrap() {
                new.insert(String::from(key), val.to_owned());
            }
        }
    }
    new
}
