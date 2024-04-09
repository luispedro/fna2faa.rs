use bio::io::fasta;
use std::io;
use std::io::Write;


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


    /// Output file (if not specified, output to stdout)
    #[arg(short, long)]
    output: Option<String>,

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
    let writer_r : Box<dyn Write> = match args.output {
        None => Box::new(io::stdout()),
        Some(s) => if s == "-" {
                Box::new(io::stdout())
            } else {
                Box::new(std::fs::File::create(s).unwrap())
            }
    };
    let mut writer = io::BufWriter::new(writer_r);

    let mut prot = Vec::new();
    let mut rseq = Vec::new();

    for record in reader.records() {
        let record = record.unwrap();
        let fseq = record.seq();
        // It is very ugly to repeat the code, but all the solutions I tried
        // had a 5-10% performance penalty.
        if args.all_frames {
            molbio::rev_compl_to(fseq, &mut rseq);
            for frame in 0..6 {
                let seq = if frame > 2 { &rseq } else { fseq };
                prot.clear();
                let mut start_ix = frame % 3;
                while (start_ix + 3) <= seq.len() {
                    let codon = &seq[start_ix..start_ix + 3];
                    let aa = *encoded_table.get(nuc_tab.encode3(codon)).unwrap();
                    prot.push(aa);
                    if args.first_stop && aa == b'*' { break; }
                    start_ix += 3;
                }
                writer.write(b">").unwrap();
                writer.write(record.id().as_bytes()).unwrap();
                match frame {
                    0 => { writer.write(b":0:0\n").unwrap(); }
                    1 => { writer.write(b":1:0\n").unwrap(); }
                    2 => { writer.write(b":2:0\n").unwrap(); }
                    3 => { writer.write(b":0:1\n").unwrap(); }
                    4 => { writer.write(b":1:1\n").unwrap(); }
                    5 => { writer.write(b":2:1\n").unwrap(); }
                    _ => { unreachable!(); }
                }
                writer.write(&prot).unwrap();
                writer.write(b"\n").unwrap();

            }
        } else {
            let seq = if rc { molbio::rev_compl_to(fseq, &mut rseq); &rseq } else { fseq };
            prot.clear();
            let mut start_ix = frame;
            while (start_ix + 3) <= seq.len() {
                let codon = &seq[start_ix..start_ix + 3];
                let aa = *encoded_table.get(nuc_tab.encode3(codon)).unwrap();
                prot.push(aa);
                if args.first_stop && aa == b'*' { break; }
                start_ix += 3;
            }
            writer.write(b">").unwrap();
            writer.write(record.id().as_bytes()).unwrap();
            writer.write(b"\n").unwrap();
            writer.write(&prot).unwrap();
            writer.write(b"\n").unwrap();
        }
    }
    writer.flush().unwrap();
}

