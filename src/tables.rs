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
    pub fn encode3(&self, cod: &[u8]) -> usize {
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

fn ambiguous_codon2aa(map: &HashMap<char, &'static str>, codon: &str) -> Option<char> {
    let mut result = None;

    if codon.len() != 3 {
        return result;
    }

    let chars: Vec<char> = codon.chars().collect();
    for n1 in map[&chars[0]].chars() {
        for n2 in map[&chars[1]].chars() {
            for n3 in map[&chars[2]].chars() {
                let nc = format!("{}{}{}", n1, n2, n3);
                if let Some(aa) = codon2aa().get(&nc[..]) {
                    if result.is_none() {
                        result = Some(*aa);
                    } else if result.unwrap() != *aa {
                        return None;
                    }
                }
            }
        }
    }

    result
}

#[test]
fn test_ambiguous2aa() {
    assert!(ambiguous_codon2aa(&ambiguous2nucs(), "MGR") == Some('R'));
    assert!(ambiguous_codon2aa(&ambiguous2nucs(), "TGN") == None);
    assert!(ambiguous_codon2aa(&ambiguous2nucs(), "TGY") == Some('C'));
    assert!(ambiguous_codon2aa(&ambiguous2nucs(), "TGT") == Some('C'));
    assert!(ambiguous_codon2aa(&ambiguous2nucs(), "TGC") == Some('C'));
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

fn general_codon2aa(ambiguous_map: &HashMap<char, &'static str>,
                    codon2aa_tab: &HashMap<&'static str, char>,
                    c : &str) -> Option<char> {
    if c.chars().all(|ch| "ACGT".contains(ch)) {
        codon2aa_tab.get(&c).map(|&aa| aa)
    } else {
        ambiguous_codon2aa(&ambiguous_map, &c)
    }
}


pub fn build_table(nuc_tab : &CodonEncoder) -> Vec<u8> {
    let mut vec = vec![b'X'; 16*16*16];
    let codon2aa_tab = codon2aa();
    let ambiguous_map = ambiguous2nucs();

    for c in generate_all_codons(ambiguous_map.keys().cloned().collect()) {
        if let Some(aa) = general_codon2aa(&ambiguous_map, &codon2aa_tab, &c) {
            vec[nuc_tab.encode3(c.as_bytes())] = aa as u8;
        }
    }
    vec
}


