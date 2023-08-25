/// Reverse complement a DNA string
pub fn rev_compl_to(seq: &[u8], res: &mut Vec<u8>) {
    res.clear();
    res.reserve(seq.len());
    for c in seq.iter().rev() {
        res.push(match c {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',

            // Uracyl becomes A
            b'U' => b'A',

            // Handle the ambiguous ones
            b'M' => b'K',
            b'R' => b'Y',
            b'W' => b'W',
            b'S' => b'S',
            b'Y' => b'R',
            b'K' => b'M',
            b'V' => b'B',
            b'H' => b'D',
            b'D' => b'H',
            b'B' => b'V',
            b'X' => b'X',
            b'N' => b'N',

            _ => *c,
        });
    }
}



#[test]
fn test_rc() {
    fn rc(seq: &[u8]) -> Vec<u8> {
        let mut res = Vec::with_capacity(seq.len());
        rev_compl_to(seq, &mut res);
        res
    }
    assert!(rc(b"ATTC") == b"GAAT");
    assert!(rc(b"TTNN") == b"NNAA");
    assert!(rc(b"AUMW") == b"WKAT");
}

