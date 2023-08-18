use std::collections::HashMap;

pub struct CodonEncoder {
    tab: Vec<usize>,
}
impl CodonEncoder {
    pub fn mk_encoder() -> CodonEncoder {
        let mut vec = vec![0; 256];
        for c in "ATCGMRWSYKVHDBXUN".chars() {
            vec[c as usize] = encode1slow(c).unwrap();
            let c = c.to_ascii_lowercase();
            vec[c as usize] = encode1slow(c).unwrap();
        }
        CodonEncoder { tab: vec }
    }
    pub fn encode3(&self, n: &str) -> usize {
        let cod = n.as_bytes();
        if cod.len() < 3 {
            0
        } else {
            16 * 16 * self.tab[usize::from(cod[0])] +
                 16 * self.tab[usize::from(cod[1])] +
                      self.tab[usize::from(cod[2])]
        }
    }
}


fn codon2aa() -> HashMap<&'static str, char> {
    let mut map = HashMap::new();
    map.insert("TTT", 'F');
    map.insert("TTC", 'F');
    map.insert("TTA", 'L');
    map.insert("TTG", 'L');

    map.insert("TCT", 'S');
    map.insert("TCC", 'S');
    map.insert("TCA", 'S');
    map.insert("TCG", 'S');

    map.insert("TAT", 'Y');
    map.insert("TAC", 'Y');
    map.insert("TAA", '*');
    map.insert("TAG", '*');

    map.insert("TGT", 'C');
    map.insert("TGC", 'C');
    map.insert("TGA", '*');
    map.insert("TGG", 'W');

    map.insert("CTT", 'L');
    map.insert("CTC", 'L');
    map.insert("CTA", 'L');
    map.insert("CTG", 'L');

    map.insert("CCT", 'P');
    map.insert("CCC", 'P');
    map.insert("CCA", 'P');
    map.insert("CCG", 'P');

    map.insert("CAT", 'H');
    map.insert("CAC", 'H');
    map.insert("CAA", 'Q');
    map.insert("CAG", 'Q');

    map.insert("CGT", 'R');
    map.insert("CGC", 'R');
    map.insert("CGA", 'R');
    map.insert("CGG", 'R');

    map.insert("ATT", 'I');
    map.insert("ATC", 'I');
    map.insert("ATA", 'I');
    map.insert("ATG", 'M');

    map.insert("ACT", 'T');
    map.insert("ACC", 'T');
    map.insert("ACA", 'T');
    map.insert("ACG", 'T');

    map.insert("AAT", 'N');
    map.insert("AAC", 'N');
    map.insert("AAA", 'K');
    map.insert("AAG", 'K');

    map.insert("AGT", 'S');
    map.insert("AGC", 'S');
    map.insert("AGA", 'R');
    map.insert("AGG", 'R');

    map.insert("GTT", 'V');
    map.insert("GTC", 'V');
    map.insert("GTA", 'V');
    map.insert("GTG", 'V');

    map.insert("GCT", 'A');
    map.insert("GCC", 'A');
    map.insert("GCA", 'A');
    map.insert("GCG", 'A');

    map.insert("GAT", 'D');
    map.insert("GAC", 'D');
    map.insert("GAA", 'E');
    map.insert("GAG", 'E');

    map.insert("GGT", 'G');
    map.insert("GGC", 'G');
    map.insert("GGA", 'G');
    map.insert("GGG", 'G');

    map
}


fn ambiguous2nucs() -> HashMap<char, &'static str> {
    let mut map = HashMap::new();
    map.insert('M', "AC");
    map.insert('R', "AG");
    map.insert('W', "AT");
    map.insert('S', "CG");
    map.insert('Y', "CT");
    map.insert('K', "GT");
    map.insert('V', "ACG");
    map.insert('H', "ACT");
    map.insert('D', "AGT");
    map.insert('B', "CGT");
    map.insert('X', "GACT");
    map.insert('N', "GACT");

    // These are convenient to simplify the code even though they are not
    // ambiguous.
    map.insert('U', "T");
    map.insert('A', "A");
    map.insert('C', "C");
    map.insert('G', "G");
    map.insert('T', "T");
    map
}

fn expand_all_ambiguous_nucs(map: &HashMap<char, &'static str>, codon: &str) -> Vec<String> {
    let mut results = Vec::new();

    if codon.len() != 3 {
        return results;
    }

    let chars: Vec<char> = codon.chars().collect();
    for n1 in map[&chars[0]].chars() {
        for n2 in map[&chars[1]].chars() {
            for n3 in map[&chars[2]].chars() {
                results.push(format!("{}{}{}", n1, n2, n3));
            }
        }
    }

    results
}


fn encode1slow(n: char) -> Option<usize> {
    let n = n.to_ascii_uppercase();
    let n = if n == 'U' { 'T' } else { n };
    let n = if n == 'N' { 'X' } else { n };
    "ATCGMRWSYKVHDBX".find(n)
}


fn generate_all_codons(alphabet : Vec<char>) -> Vec<String> {
    let mut codons = Vec::new();

    for n1 in &alphabet {
        for n2 in &alphabet {
            for n3 in &alphabet {
                codons.push(format!("{}{}{}", n1, n2, n3));
            }
        }
    }

    codons
}


pub fn build_table(nuc_tab : &CodonEncoder) -> Vec<char> {
    let mut vec = vec!['X'; 16*16*16];
    let codon2aa_fn = codon2aa();
    let ambiguous_map = ambiguous2nucs();

    for c in generate_all_codons(ambiguous_map.keys().cloned().collect()) {
        if c.chars().all(|ch| "ACGT".contains(ch)) {
            if let Some(&aa) = codon2aa_fn.get(&c.as_str()) {
                vec[nuc_tab.encode3(&c)] = aa;
            }
        } else {
            let mut prev: Option<char> = None;
            let mut all_same = true;

            for n in expand_all_ambiguous_nucs(&ambiguous_map, &c) {
                if let Some(&aa) = codon2aa_fn.get(&n.as_str()) {
                    if let Some(prev_a) = prev {
                        if aa != prev_a {
                            all_same = false;
                            break;
                        }
                    } else {
                        prev = Some(aa);
                    }
                }
            }

            if all_same {
                if let Some(aa) = prev {
                    vec[nuc_tab.encode3(&c)] = aa;
                }
            }
        }
    }
    vec
}


