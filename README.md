# make_prg

A tool to create a PRG for input to [Pandora][pandora] and [Gramtools][gramtools] from a
Multiple Sequence Alignment.

[TOC]: #

## Table of Contents
- [Install](#install)
  - [Conda](#conda)
  - [Local](#local)
  - [Container](#container)
- [Usage](#usage)
  - [CLI](#cli)
  - [Nextflow](#nextflow)
- [Input](#input)
- [Changing parameters](#changing-parameters)

## Install

### Conda

[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/make_prg)](https://anaconda.org/bioconda/make_prg)
[![bioconda version](https://anaconda.org/bioconda/make_prg/badges/platforms.svg)](https://anaconda.org/bioconda/make_prg)

Prerequisite: [`conda`][conda] (and bioconda channel [correctly set up][channels])

```sh
conda install make_prg
```

### Local

Requirements: `python>=3`

```sh
git clone https://github.com/rmcolq/make_prg.git
cd make_prg
python -m pip install .
make_prg --help
```

To additionally run the tests

```shell
python -m pip install nose hypothesis
nosetests tests/
```

This installs the CLI tool `make_prg`. The nextflow script `make_prg_nexflow.nf` assumes
that `make_prg` is installed.

### Container

[![Docker Repository on Quay](https://quay.io/repository/iqballab/make_prg/status "Docker Repository on Quay")](https://quay.io/repository/iqballab/make_prg)

Containers for this tool are [hosted on quay.io][tags].

An example, of running it in a [Singularity][singularity] container would be

```
tag="latest"
URI="docker://quay.io/iqballab/make_prg:${tag}"
singularity exec "$URI" make_prg --help
```

A list of the valid tags can be found [here][tags].

## Usage

### CLI

```
$ make_prg --help
usage: make_prg <subcommand> <options>

Subcommand entrypoint

optional arguments:
  -h, --help     show this help message and exit
  -V, --version  show program's version number and exit
  -v, --verbose  Run with high verbosity (debug level logging)

Available subcommands:

    from_msa     Make PRG from multiple sequence alignment
```

#### `from_msa`

```
$ make_prg from_msa --help
usage: make_prg from_msa [options] <MSA input file>

positional arguments:
  MSA                   Input file: a multiple sequence alignment

optional arguments:
  -h, --help            show this help message and exit
  -f ALIGNMENT_FORMAT, --alignment_format ALIGNMENT_FORMAT
                        Alignment format of MSA, must be a biopython AlignIO input alignment_format. See http://biopython.org/wiki/AlignIO. Default: fasta
  -N MAX_NESTING, --max_nesting MAX_NESTING
                        Maximum number of levels to use for nesting. Default: 5
  -L MIN_MATCH_LENGTH, --min_match_length MIN_MATCH_LENGTH
                        Minimum number of consecutive characters which must be identical for a match. Default: 7
  -o OUTPUT_DIR, --outdir OUTPUT_DIR
                        Output directory. Default: .
  -n PRG_NAME, --prg_name PRG_NAME
                        Prg file name. Default: MSA file name
  -S SEQID, --seqid SEQID
                        Sequence identifier to use for the output sequence/PRG. Default is the file name
  --no_overwrite        Do not replace an existing prg file
  --summary             Write a summary file
  -O OUTPUT_TYPE, --output-type OUTPUT_TYPE
                        p: PRG, b: Binary, g: GFA, a: All. Combinations are allowed i.e., gb: GFA and Binary. Default: a
  --log LOG             Path to write log to. Default is stderr
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

[channels]: https://bioconda.github.io/user/install.html#set-up-channels
[gramtools]: https://github.com/iqbal-lab-org/gramtools
[nf]: https://www.nextflow.io/
[pandora]: https://github.com/rmcolq/pandora
[singularity]: https://sylabs.io/
[tags]: https://quay.io/repository/iqballab/make_prg?tab=tags
[conda]: https://conda.io
