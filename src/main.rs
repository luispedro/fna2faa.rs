pub mod tables;

fn main() {
    let nuc_tab = tables::CodonEncoder::mk_encoder();
    let encoded_table = tables::build_table(&nuc_tab);
    for aa in encoded_table {
        print!("{}", aa);
    }
    println!();
}

