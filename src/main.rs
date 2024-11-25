use bio::io::fasta;
use std::env;
fn main() -> () {
    let mut arguments= Vec::new();
    for argument in env::args() {
        arguments.push(argument); 
    }
    let path_to_file= &arguments[1];
    let sh = match fasta::Reader::from_file(path_to_file) {
        Ok(val) => val,
        Err(_e) => panic!("could not read fasta file"),
    };
    let mut gcount = 0i64;
    let mut ccount = 0i64;
    let mut allcount = 0i64;
    for result in sh.records() {
        let record = result.unwrap();
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
        let gcprop = (gcount + ccount) as f64 / allcount as f64;
        println!(
            "Number of Gs\tNumber of Cs\tTotal number\tProportion of GC\n{gcount}\t{ccount}\t{allcount}\t{gcprop}"
        )

    }
}
