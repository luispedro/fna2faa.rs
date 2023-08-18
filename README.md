# Convert nucleotide FASTA to amino acid FASTA

Very fast conversion of nucleotide to amino acids FASTA. Correctly handles ambiguous nucleotides in IUPAC notation (_i.e._, if if all possible underlying nucleotides translate to the same codon, then this codon is output).

## Speed

Compared to the [C++ implementation](https://github.com/luispedro/fna2faa), this is >10x faster. It is even 7x than the C++ version using the same approach that does not support ambiguous nucleotides.

| Test       | Rust (this one)  | C++ (full) | C++ (no ambiguous) |
| ---------- | ---------------- | ---------- | ------------------ |
| 2.5M genes | 0.85s            | 11s32      |    6s37            |


- _Author_: [Luis Pedro Coelho](http://luispedro.org)
- _License_: MIT

