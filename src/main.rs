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

    /// Stop translation when first stop codon is found
    #[arg(short = 's', long)]
    first_stop: bool,

    /// Output all frames
    #[arg(short, long)]
    all_frames : bool,

}



fn main() {
    let args = Args::parse();
    let fname = if args.fname == "-" { "/dev/stdin" } else { &args.fname };

    args.frame.and_then(|f| if f > 5 { panic!("Invalid frame (should be 0..5)") } else { Some(()) });
    let frame = args.frame.unwrap_or(0) % 3;
    let rc = args.frame.unwrap_or(0) > 2;

    let nuc_tab = tables::CodonEncoder::mk_encoder();
    let encoded_table = tables::build_table(&nuc_tab);

    let reader = fasta::Reader::new(io::BufReader::new(std::fs::File::open(fname).unwrap()));
    let mut writer = fasta::Writer::new(io::stdout());

    let mut prot = Vec::new();
    let mut rseq = Vec::new();

    for record in reader.records() {
        let record = record.unwrap();
        let dseq = record.seq().to_vec();
        // It is very ugly to repeat the code, but all the solutions I tried
        // had a 5-10% performance penalty.
        if args.all_frames {
            molbio::rev_compl_to(record.seq(), &mut rseq);
            for frame in 0..6 {
                let seq = if frame > 2 { &rseq } else { &dseq };
                prot.clear();
                let mut start_ix = frame % 3;
                while (start_ix + 3) <= seq.len() {
                    let codon = &seq[start_ix..start_ix + 3];
                    let aa = *encoded_table.get(nuc_tab.encode3(codon)).unwrap();
                    prot.push(aa);
                    if args.first_stop && aa == b'*' { break; }
                    start_ix += 3;
                }
                let nid = format!("{}:{}:{}", record.id(), frame % 3, if frame > 2 { "1" } else { "0" });
                let nrecord = fasta::Record::with_attrs(&nid, None, &prot[..]);
                writer.write_record(&nrecord).unwrap();
            }
        } else {
            let seq = if rc { molbio::rev_compl_to(record.seq(), &mut rseq); &rseq } else { &dseq };
            prot.clear();
            let mut start_ix = frame;
            while (start_ix + 3) <= seq.len() {
                let codon = &seq[start_ix..start_ix + 3];
                let aa = *encoded_table.get(nuc_tab.encode3(codon)).unwrap();
                prot.push(aa);
                if args.first_stop && aa == b'*' { break; }
                start_ix += 3;
            }
            let nrecord = fasta::Record::with_attrs(record.id(), None, &prot[..]);
            writer.write_record(&nrecord).unwrap();
        }
    }
}

