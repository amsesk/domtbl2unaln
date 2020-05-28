use crate::DOMTBL_COLUMNS;

use csv::WriterBuilder;
use multimap::MultiMap;
use std::fmt;
use std::path::Path;

#[derive(Clone)]
pub struct Hits {
    pub file_path: String,
    pub n_hits: usize,
    pub inner: MultiMap<String, Hit>,
}
impl Hits {
    pub fn new(file_path: String, inner: MultiMap<String, Hit>) -> Hits {
        Hits {
            file_path,
            n_hits: inner.len(),
            inner,
        }
    }
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
    pub fn best_filter(&mut self) {
        for (_marker, hitlist) in self.inner.iter_all_mut() {
            hitlist.sort_by(|a, b| b.fs_score.partial_cmp(&a.fs_score).unwrap());
            hitlist.truncate(1);
        }
    }
    pub fn hit_report_csv(&self, handle: &std::fs::File) -> Result<(), std::io::Error> {
        // Write all hit information to csv
        let mut wtr = WriterBuilder::new().delimiter(b'\t').from_writer(handle);
        for (_marker, hits) in self.inner.iter_all() {
            for h in hits {
                assert_eq!(_marker, &h.query);
                wtr.write_record(&[
                    &format!("{}", h.target),
                    &format!("{}", h.query),
                    &format!("{:e}", h.fs_evalue),
                    &format!("{}", h.fs_score),
                ])?;
            }
        }
        wtr.flush()?;
        Ok(())
    }
    pub fn recovery_csv(&self, handle: &std::fs::File) -> Result<(), std::io::Error> {
        let mut wtr = WriterBuilder::new().delimiter(b'\t').from_writer(handle);
        wtr.write_record(&[
            &format!("{}", self.file_path),
            &format!("{}", self.n_hits),
            &format!("{:.2}", self.perc_duplicated()),
        ])?;
        Ok(())
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
        write!(f, "{}-{:e}", self.target, self.fs_evalue)
    }
}
