# make_prg

A tool to create and update PRGs for input to [Pandora][pandora] and [Gramtools][gramtools] from a set of 
Multiple Sequence Alignments.

# Support

We fully support `make_prg` on `linux` with `python` versions `3.8`-`3.11`.
No guarantees are made for other operating systems or `python` versions.

[TOC]: #

## Table of Contents
- [Dependencies](#dependencies)
- [Install](#install)
  - [No installation needed - precompiled portable binary](#no-installation-needed---precompiled-portable-binary)
  - [pip](#pip)
  - [conda](#conda)
- [Running on a sample example](#running-on-a-sample-example)
- [Usage](#usage)

## Dependencies

`make_prg` has two commands: `from_msa` and `update`. The `update` command requires `MAFFT` to be installed:
  1. from source: https://mafft.cbrc.jp/alignment/software/;
  2. using `conda`: `conda install -c bioconda mafft`;

## Install

### No installation needed - precompiled portable binary

> Note: We provide precompiled binaries for Linux OS only

You can use `make_prg` with no installation at all by simply downloading the precompiled binary, and running it.
In this binary, all libraries are linked statically. Compilation is done using [PyInstaller](https://github.com/pyinstaller/pyinstaller).

#### Requirements
`GLIBC >= 2.27` (present on `Ubuntu >= 18.04`, `Debian >= 10`, `CentOS >= 8`, etc);

#### Download
```
wget https://github.com/iqbal-lab-org/make_prg/releases/download/0.4.0/make_prg_0.4.0
```

#### Run
```
chmod +x make_prg_0.4.0
./make_prg_0.4.0 -h
```

### pip

**Requirements**: `python>=3.8,<=3.11`

```sh
pip install make_prg
```

### conda

```sh
conda install -c bioconda make_prg
```

## Running on a sample example

To see how to input files to both `make_prg from_msa` and `make_prg update`, and the outputs
they create on a sample example, see [sample example](sample_example).

## Usage

```
$ make_prg --help
usage: make_prg <subcommand> <options>

Subcommand entrypoint

optional arguments:
  -h, --help     show this help message and exit
  -V, --version  show program's version number and exit

Available subcommands:
  
    from_msa     Make PRG from multiple sequence alignment
    update       Update PRGs given new sequences.
```

#### `from_msa`

```
$ make_prg from_msa --help
usage: make_prg from_msa

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Multiple sequence alignment file or a directory containing such files
  -s SUFFIX, --suffix SUFFIX
                        If the input parameter (-i, --input) is a directory, then filter for files with this suffix. If this parameter is not given, all files in the
                        input directory is considered.
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix for the output files
  -f ALIGNMENT_FORMAT, --alignment-format ALIGNMENT_FORMAT
                        Alignment format of MSA, must be a biopython AlignIO input alignment_format. See http://biopython.org/wiki/AlignIO. Default: fasta
  -N MAX_NESTING, --max-nesting MAX_NESTING
                        Maximum number of levels to use for nesting. Default: 5
  -L MIN_MATCH_LENGTH, --min-match-length MIN_MATCH_LENGTH
                        Minimum number of consecutive characters which must be identical for a match. Default: 7
  -O OUTPUT_TYPE, --output-type OUTPUT_TYPE
                        p: PRG, b: Binary, g: GFA, a: All. Combinations are allowed i.e., gb: GFA and Binary. Default: a
  -F, --force           Force overwrite previous output
  -t THREADS, --threads THREADS
                        Number of threads. 0 will use all available. Default: 1
  -v, --verbose         Increase output verbosity (-v for debug, -vv for trace - trace is for developers only)
  --log LOG             Path to write log to. Default is stderr
```

#### `update`

```
$ make_prg update --help
usage: make_prg update

optional arguments:
  -h, --help            show this help message and exit
  -u UPDATE_DS, --update-DS UPDATE_DS
                        Filepath to the update data structures (a *.update_DS.zip file created from make_prg from_msa or update)
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix for the output files
  -d DENOVO_PATHS, --denovo-paths DENOVO_PATHS
                        Filepath containing denovo sequences. Should point to a denovo_paths.txt file
  -D LONG_DELETION_THRESHOLD, --deletion-threshold LONG_DELETION_THRESHOLD
                        Ignores long deletions of the given size or longer. If long deletions should not be ignored, put a large value. Default: 10
  -m MAFFT, --mafft MAFFT
                        Path to MAFFT executable. By default, it is assumed to be on $PATH
  -O OUTPUT_TYPE, --output-type OUTPUT_TYPE
                        p: PRG, b: Binary, g: GFA, a: All. Combinations are allowed i.e., gb: GFA and Binary. Default: a
  -F, --force           Force overwrite previous output
  -t THREADS, --threads THREADS
                        Number of threads. 0 will use all available. Default: 1
  -v, --verbose         Increase output verbosity (-v for debug, -vv for trace - trace is for developers only)
  --log LOG             Path to write log to. Default is stderr
```

[pandora]: https://github.com/rmcolq/pandora
[gramtools]: https://github.com/iqbal-lab-org/gramtools/
