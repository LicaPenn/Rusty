use bio::io::fasta;

fn main() -> () {



    let sh = match fasta::Reader::from_file("./data/sequence.fasta") {
        Ok(val) => val,
        Err(_e) => panic!("could not read fasta file"),
    };
    let mut gcount = 0;
    let mut ccount= 0;
    let mut allcount=0; 
    for result in sh.records() {
        let record = result.unwrap();
        let sequence = record.seq();
        for letter in sequence {
            if letter == &71u8 {gcount +=1} 
            if letter == &67u8 {ccount +=1} 
            allcount +=1;

        }
   let gcprop= (gcount+ccount) as f64 / allcount as f64;
        println!("Number of Gs: {}, number of Cs: {}, total number:{}, proportion of gc:{}", gcount, ccount, allcount,gcprop)
    }

 
}

