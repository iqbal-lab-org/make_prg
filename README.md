# make_prg

A tool to create a PRG for input to [Pandora][pandora] from a Multiple Sequence
Alignment.

[TOC]: #

# Table of Contents
- [Install](#install)
  - [Local](#local)
  - [Singularity](#singularity)
- [Usage](#usage)
  - [CLI](#cli)
  - [Nextflow](#nextflow)
- [Input](#input)
- [Changing parameters](#changing-parameters)


## Install

### Local

Requirements: `python>=3`

```sh
git clone https://github.com/rmcolq/make_prg.git
cd make_prg
python3 setup.py test
python3 setup.py install
```

This installs the CLI tool `make_prg`. The nextflow script `make_prg_nexflow.nf` assumes
that `make_prg` is installed.

### Singularity

Requirement: `singularity>=3`

```sh
URI="library://mbhall88/default/make_prg"
singularity exec "$URI" make_prg --help
```

You can also get a build for a specific commit. A list of the valid tags can be found
[here][tags].

```sh
tag="0bb4a27"
URI="library://mbhall88/default/make_prg:$tag"
singularity exec "$URI" make_prg --help
```

## Usage

### CLI

```
$ make_prg --help
usage: make_prg <subcommand> <options>

script to run make_prg subcommands

optional arguments:
  -h, --help    show this help message and exit
  --version     show program's version number and exit

Available subcommands:

    prg_from_msa
                Make PRG from multiple sequence alignment
```

#### `prg_from_msa`

```
$ make_prg prg_from_msa --help
usage: make_prg prg_from_msa [options] <MSA input file>

positional arguments:
  MSA                   Input file: a multiple sequence alignment in supported alignment_format. If not in aligned fasta
                        alignment_format, use -f to input the alignment_format type

optional arguments:
  -h, --help            show this help message and exit
  -f ALIGNMENT_FORMAT, --alignment_format ALIGNMENT_FORMAT
                        Alignment format of MSA, must be a biopython AlignIO input alignment_format. See
                        http://biopython.org/wiki/AlignIO. Default: fasta
  --max_nesting MAX_NESTING
                        Maximum number of levels to use for nesting. Default: 5
  --min_match_length MIN_MATCH_LENGTH
                        Minimum number of consecutive characters which must be identical for a match. Default: 7
  -p OUTPUT_PREFIX, --prefix OUTPUT_PREFIX
                        Output prefix
  --no_overwrite        Do not overwrite pre-existing prg file with same name
  -v, --verbose         Run with high verbosity (debug level logging)
```

### Nextflow

Requirements: [Nextflow][nf]

```
    Usage: nextflow run make_prg_nexflow.nf <arguments>

    Required arguments:
      --tsv_in  FILENAME  An index file of MSA to build PRGs of

    Optional arguments:
```

## Input

Multiple Sequence Alignment files for genes/dna sequences for which we wantPRGs, and an
tab-separated index of these in the form:

```
sample_id       infile
GC0000001   /absolute/path/to/GC0000001_na_aln.fa.gz
GC0000002   /absolute/path/to/GC0000002_na_aln.fa
```

## Changing parameters

There are some parameters at the top of the nextflow file which could be changed but
which I have not made command line parameters:

```
max_nesting             This is the maximum number depth of bubbles in PRG, setting to 1 will allow variants, \\
                        but no nesting
min_match_length        Controls graph complexity
alignment_format        Any format accepted by biopython's AlignIO
max_forks_make_prg      If working on a cluster which allows unlimited parallel jobs per user, this will be \\
                        used by nextflow to control maximum number of processes of this type that can run in \\
                        parallel.
max_forks_make_fasta
```

[tags]: https://cloud.sylabs.io/library/mbhall88/default/make_prg
[pandora]: https://github.com/rmcolq/pandora
[nf]: https://www.nextflow.io/
