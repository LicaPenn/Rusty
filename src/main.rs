use anyhow::{Context, Result};
use bio::io::fasta::{self, Reader};
use clap::{arg, command, value_parser};
use std::{fs::File, io::BufReader, path::PathBuf};

fn main() -> Result<()> {
    let matches = command!()
        .about("Counts G and C of your fasta file")
        .arg(
            arg!(
                 <FILE> "Path to Fasta file"
            )
            .required(true)
            .value_parser(value_parser!(PathBuf)),
        )
        .get_matches();

    let path_to_file = matches.get_one::<PathBuf>("FILE").unwrap();

    let sh =
        fasta::Reader::from_file(path_to_file).context("Did you spell the path name correctly?")?;

    let (gcount, ccount, allcount, gcprop) = gccalculator(sh)?;

    println!(
            "Number of Gs\tNumber of Cs\tTotal number\tProportion of GC\n{gcount}\t{ccount}\t{allcount}\t{gcprop}"
        );

    Ok(())
}

fn gccalculator(sh: Reader<BufReader<File>>) -> Result<(i64, i64, i64, f64)> {
    let mut gcount = 0i64;
    let mut ccount = 0i64;
    let mut allcount = 0i64;
    for result in sh.records() {
        let record = result.context("Could not read record.")?;
        let sequence = record.seq();
        for letter in sequence {
            if letter == &71u8 {
                gcount += 1;
            }
            if letter == &67u8 {
                ccount += 1;
            }
            allcount += 1;
        }
    }
    let gcprop = (gcount + ccount) as f64 / allcount as f64;
    Ok((gcount, ccount, allcount, gcprop))
}

#[cfg(test)]
mod tests {

    const DATA: &str = "GGGCCCGGGCC";

    #[test]
    fn test_gc() {
        // logic
    }
}
