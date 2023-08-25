use bio::io::fasta;
use std::io;

pub mod tables;
pub mod molbio;

use clap::Parser;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg()]
    fname: String,

    /// Specify which frame to output (0, 1, 2, 3, 4, 5)
    /// Frames 0-2 are -> strand, 3-5 are <-.
    #[arg(short, long)]
    frame: Option<usize>,

}

fn main() {
    let args = Args::parse();
    let fname = if args.fname == "-" { "/dev/stdin" } else { &args.fname };

    args.frame.and_then(|f| if f > 5 { panic!("Invalid frame (should be 0..5)") } else { Some(()) });
    let frame = args.frame.unwrap_or(0) % 3;
    let use_complement = args.frame.unwrap_or(0) >= 3;

    let nuc_tab = tables::CodonEncoder::mk_encoder();
    let encoded_table = tables::build_table(&nuc_tab);

    let reader = fasta::Reader::new(io::BufReader::new(std::fs::File::open(fname).unwrap()));
    let mut writer = fasta::Writer::new(io::stdout());

    let mut prot = Vec::new();
    for record in reader.records() {
        prot.clear();
        let record = record.unwrap();
        let seq = if use_complement { molbio::rev_compl(record.seq()) } else { record.seq().to_vec() };
        let mut start_ix = frame;
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

