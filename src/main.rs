use anyhow::{Context, Result};
use clap::{arg, command, value_parser};
use needletail::{parse_fastx_file, FastxReader, Sequence};
use std::path::PathBuf;

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

    let sh = parse_fastx_file(path_to_file).context("Did you spell the path name correctly?")?;

    let (gcount, ccount, allcount, gcprop) = gccalculator(sh)?;

    println!(
            "Number of Gs\tNumber of Cs\tTotal number\tProportion of GC\n{gcount}\t{ccount}\t{allcount}\t{gcprop}"
        );

    Ok(())
}

fn gccalculator(mut sh: Box<dyn FastxReader>) -> Result<(i64, i64, i64, f64)> {
    let mut gcount = 0i64;
    let mut ccount = 0i64;
    let mut allcount = 0i64;
    while let Some(record) = sh.next() {
        let record = record?;
        let sequence = record.sequence();
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
    use super::*;

    const DATA: &str = ">test\nGGGCCCGGGCCC";

    #[test]
    fn test_gc() {
        let reader = needletail::parse_fastx_reader(DATA.as_bytes()).unwrap();
        let (gcount, ccount, allcount, gcprop) = gccalculator(reader).unwrap();

        assert_eq!(gcount, 6);
        assert_eq!(ccount, 6);
        assert_eq!(allcount, 12);
        assert_eq!(gcprop, 1.0);
    }
}
