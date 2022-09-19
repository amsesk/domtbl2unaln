extern crate bio;
use bio::utils::Text;
use clap::{App, Arg};
use domtbl2unaln::hits::Hits;
use domtbl2unaln::{file_list, has_enough_occupants, parse_and_filter, parse_cutoffs};
use std::boxed::Box;
use std::fs::OpenOptions;

use std::fs::{read_to_string, File};

use std::io::BufWriter;
use std::io::Write;
use std::path::Path;

pub fn main() -> Result<(), std::io::Error> {
    let args = App::new("domtbl-reader")
        .version("0.1")
        .author("Kevin Amses")
        .about("Reads and parses domtbl from hmmsearch.")
        .arg(
            Arg::new("domtbls")
                .long("domtbls")
                .value_name("DIR")
                .about("Path to directory containing domtbls.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::new("proteins")
                .long("proteins")
                .value_name("DIR")
                .about("A protein fasta including predictions for ALL taxa, with samtools faidx.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::new("outdir")
                .long("outdir")
                .value_name("DIR")
                .about("Path to existing directory to dump unaligned fastas.")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::new("best")
                .long("best")
                .about("Pass this flag to take the best hits for each marker only.")
                .takes_value(false)
                .required(false),
        )
        .arg(
            Arg::new("cutoffs")
                .long("cutoffs")
                .about("Path to busco cutoffs file.")
                .takes_value(true)
                .required(false),
        )
        .arg(
            Arg::new("occupancy_cutoff")
                .long("occupancy_cutoff")
                .about("Marker occupancy cutoff for not printing unaln. Default = 0.75")
                .default_value("0.75")
                .takes_value(true)
                .required(false),
        )
        .get_matches();

    // Unwrap on None here is impossible because
    // clap will have already killed the process
    // since a required argument is missing
    let domtbls = Path::new(args.value_of("domtbls").unwrap());
    let proteins = Path::new(args.value_of("proteins").unwrap())
        .canonicalize()
        .unwrap();
    let outdir = Path::new(args.value_of("outdir").unwrap())
        .canonicalize()
        .unwrap();
    let occ_cutoff: f64 = args.value_of("occupancy_cutoff").unwrap().parse().unwrap();

    let busco_filt;
    let mut cutoffs = domtbl2unaln::FUNGI_ODB10_CUTOFFS;
    match args.value_of("cutoffs") {
        Some(c) => {
            match c {
                "fungi" => (),
                "mollicutes" => cutoffs = domtbl2unaln::MOLLICUTES_ODB10_CUTOFFS,
                "burkholderiales" => cutoffs = domtbl2unaln::BURKHOLDERIALES_ODB10_CUTOFFS,
                "jgi434" => cutoffs = domtbl2unaln::JGI434_ODB10_CUTOFFS,
                "mucorales" => cutoffs = domtbl2unaln::MUCORALES_ODB10_CUTOFFS,
                "bacteria" => cutoffs = domtbl2unaln::BACTERIA_ODB10_CUTOFFS,
                _ => {
                    let mpath = Path::new(&c).canonicalize().unwrap();
                    match mpath.exists() {
                        true => {
                            //let dest = "../lib/custom_cutoffs";
                            //println!("{:?} | {:?}", &mpath, &dest);
                            //copy(&mpath, &dest).unwrap();
                            //pub static CUSTOM_CUTOFFS: &str = include_str!("../lib/custom_cutoffs");
                            //cutoffs = CUSTOM_CUTOFFS
                            let cust = read_to_string(&mpath).unwrap();
                            cutoffs = Box::leak(cust.into_boxed_str());
                        }
                        false => panic!("You selected cutoffs that don't exist."),
                    };
                }
            }
            println!("Filtering by BUSCO odb10 cutoffs");
            busco_filt = true;
        }
        None => busco_filt = false,
    };

    let paths = file_list(&domtbls, "domtbl").unwrap();
    let mut all_hits: Vec<Hits> = paths
        .into_iter()
        .map(|p| parse_and_filter(p, busco_filt, cutoffs, &outdir).unwrap())
        .collect();

    // Write reports to csv
    // Truncate previous files with the same name
    let report_handle = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(outdir.join(Path::new("hit_report_all.tsv")))?;
    let recovery_handle = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(outdir.join(Path::new("marker_recovery.tsv")))?;
    for hits in all_hits.iter() {
        hits.hit_report_csv(&report_handle)?;
        hits.recovery_csv(&recovery_handle)?;
    }

    if args.occurrences_of("best") != 0 {
        println!("Doing best");
        for h in all_hits.iter_mut() {
            h.best_filter();
        }
        //hits.iter_mut().map(|(_p, h)| h.best_filter());
    }
    let odb10_marker_cutoffs = parse_cutoffs(cutoffs);

    for key in odb10_marker_cutoffs.keys() {
        println!("{}", key);
    }

    let mut index_reader = bio::io::fasta::IndexedReader::from_file(&proteins).unwrap();
    println!(
        "There are {} domtbls (and species) in this dataset.",
        all_hits.len()
    );
    for m in odb10_marker_cutoffs.keys() {
        if !has_enough_occupants(&m, &all_hits, occ_cutoff, all_hits.len() as f64) {
            continue;
        }
        let mut marker_fname = String::from(*m);
        marker_fname.push_str(".fasta");
        let outpath = outdir.join(Path::new(&marker_fname));
        let unaln = File::create(&outpath)?;
        let mut unaln_writer = BufWriter::new(&unaln);
        for hits in &all_hits {
            match hits.hits_for(&m) {
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
