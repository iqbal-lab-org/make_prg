# Toy example

Here we present a walkthrough of running `make_prg` on a toy example.
We run:
1) `make_prg from_msa` to create PRGs from MSAs;
2) `make_prg update` to update PRGs with denovo paths;

# Dependencies

* `mafft` in `PATH`. Can be installed:
  1. from source: https://mafft.cbrc.jp/alignment/software/;
  2. using `conda`: `conda install -c bioconda mafft`;

There is no need to have `make_prg` itself installed.
The test script downloads the precompiled binary of the latest release.

## Input data description

```
msas/ : contains the MSAs of 2 genes we are using as toy example here;
denovo_paths/denovo_paths.txt : contains some denovo paths on these 2 genes. This file is produced by pandora;
```

## Running

```
./run_make_prg_on_sample_example.sh
```

### Quick look at the output

### Output files from `make_prg from_msa`:
1. `msas_output/sample.prg.fa`: contains the PRGs built from the MSAs in fasta format;
2. `msas_output/sample.prg.bin.zip`: contains the PRGs built from the MSAs in a binary format;
3. `msas_output/sample.prg.gfa.zip`: contains the GFA representation of the PRGs;
4. `msas_output/sample.update_DS.zip`: contains data structures that make the PRGs updatable;

### Output files from `make_prg update`:

The output files from `make_prg update` are identical to `make_prg from_msa`, except that the PRGs are updated with
respect to the variants described in `denovo_paths/denovo_paths.txt`

### Checking the denovo variants added to the graph

Diffing `msas_output/sample.prg.fa` and `msas_updated/updated_sample.prg.fa`, we can see that the updated PRGs have more sites/alleles:

```
$ diff msas_output/sample.prg.fa msas_updated/updated_sample.prg.fa 
2c2
< TTGAGTAAAACAATCCCCCGCGCTTATATAAGCGCGTTGATATTTTTAATTATTAACAAGCAACATCATGCTAATACAGACATACAAGGAGATCATCTCTCTTTGCCTGTTTTTTATTATTTCAGGAGTGTAAACACATTTTCCG 5 T 6 C 5 CTCCCTGGCTAAT 7 C 8 A 7 ACCACATTGGCATTTATGGAGCACATCACAATATTTCAATACCATTAAAGCACTGCA 9 C 10 T 9 CAAAATGAAACACTGCGA 11 C 12 T 11 ATTAAAATT 13 A 14 C 13 TTTCAATT
---
> TTGAGTAAAACAATCCCCCGCGCTTATATAAGCGCGTTGATATTTTTA 5 A 6 G 5 TTATTAACAAGCAACATCATGCTAATACAGACATACAAGGAGATCATCTCTCTTTGCCTGTTTTTTATTATTTCAGGAGTGTAAACACATTTTCCG 7 T 8 C 7 CTCCCTGGCTAAT 9 C 10 A 9 ACCACATTGGCATTTATGGAGCACATCACAATATTTCAATACCATTAAAGCACTGCA 11 C 12 T 11 CAAAATGAAACACTGCGA 13 C 14 T 13 ATTAAAATT 15 A 16 C 15 TTTCAATT
4c4
< ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCACCGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC 5 C 6 G 5 GAAATGCCCACGACCGGGCTGGATGAGCTGACAGAGGCCGAACGCGAGAC 7 A 8 G 7 CTGGCCGAATGGCTGGCCAAGCGCCGGGAAGCCTCGCAGAAGTCGCAGGAGGCCTACACGGCCATGTCTGCGGATCGGTGGCTGGTCACGCTGGCCAAGGCCATCAGGGAAGGGCAGGA 9 GCTA 10 ACTG 9 CGCCCCGAACAGGCGGCCGCGATCTGGCACGGCATGGGGGA 11 A 12 G 11 GTCGGCAAGGCCTTGCGCAAGGCTGGTCACGCGAAGCCCAAGGCGGTCAGAAAGGGCAAGCCGGTCGATCCGGCTGATCCCAAGGATCAAGGGGAGGGGGCACCAAAGGGGAAATGA
---
> ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCA 5 C 6 T 5 CGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC 7 C 8 G 7 GAAATGCCCACGACCGGGCTGGATGAGCTGACAGAGGCCGAACGCGAGAC 9 A 10 G 9 CTGGCCGAATGGCTGGCCAAGCGCCGGGAAGCCTCGCAGAAGTCGCAGGAGGCCTACACGGCCATGTCTGCGGATCGGTGGCTGGTCACGCTGGCCAAGGCCATCAGGGAAGGGCAGGA 11 GCTA 12 ACTG 11 CGCCCCGAACAGGCGGCCGCGATCTGGCACGGCATGGGGGA 13 A 14 G 13 GTCGGCAAGGCCTTGCGCAAGGCTGGTCACGCGAAGCCCAAGGCGGTCAGAAAGGGCAAGCCGGTCGATCCGGCTGATCCCAAGGATCAAGGGG 15 A 16 T 15 GGGGGCACCAAAGGGGAAATGA
```
