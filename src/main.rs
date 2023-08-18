use bio::io::fasta;
use std::io;

pub mod tables;

fn main() {
    let nuc_tab = tables::CodonEncoder::mk_encoder();
    let encoded_table = tables::build_table(&nuc_tab);
    let reader = fasta::Reader::new(io::stdin());
    let mut writer = fasta::Writer::new(io::stdout());
    let mut prot = Vec::new();
    for record in reader.records() {
        prot.clear();
        let record = record.unwrap();
        let seq = record.seq().to_vec();
        let mut start_ix = 0;
        while (start_ix + 3) <= seq.len() {
            let codon = &seq[start_ix..start_ix + 3];
            let aa = encoded_table.get(nuc_tab.encode3(codon)).unwrap();
            prot.push(*aa);
            start_ix += 3;
        }
        let nrecord = fasta::Record::with_attrs(record.id(), None, &prot[..]);
        writer.write_record(&nrecord).unwrap();
    }
}

