extern crate bio;
use bio::utils::Text;
use clap::{App, Arg};
use domtbl2unaln::hits::Hits;
use domtbl2unaln::{file_list, has_enough_occupants, parse_and_filter, parse_cutoffs};
use std::fs::OpenOptions;

use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
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
        .arg(
            Arg::with_name("best")
                .long("best")
                .about("Pass this flag to take the best hits for each marker only.")
                .takes_value(false)
                .required(false),
        )
        .get_matches();

    // Unwrap on None here is impossible because
    // clap will have already killed the process
    // since a required argument is missing
    let domtbls = Path::new(args.value_of("domtbls").unwrap());
    let proteins = Path::new(args.value_of("proteins").unwrap());
    let outdir = Path::new(args.value_of("outdir").unwrap())
        .canonicalize()
        .unwrap();

    let paths = file_list(&domtbls, "domtbl").unwrap();
    let mut all_hits: Vec<Hits> = paths
        .into_iter()
        .map(|p| parse_and_filter(p).unwrap())
        .collect();

    // Write reports to csv
    // Truncate previous files with the same name
    let report_handle = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open("/home/aimzez/DATA/phylogeny/hit_report_all.csv")?;
    let recovery_handle = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open("/home/aimzez/DATA/phylogeny/marker_recovery.csv")?;
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

    let odb10_marker_cutoffs = parse_cutoffs(domtbl2unaln::ODB10_CUTOFFS);

    let mut index_reader = bio::io::fasta::IndexedReader::from_file(&proteins).unwrap();
    println!("{}", all_hits.len());
    for m in odb10_marker_cutoffs.keys() {
        if !has_enough_occupants(&m, &all_hits, 0.75, all_hits.len() as f64) {
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
