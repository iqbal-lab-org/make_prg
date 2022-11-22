# make_prg

A tool to create and update PRGs for input to [Pandora][pandora] from a set of 
Multiple Sequence Alignments.

[TOC]: #

## Table of Contents
- [Dependencies](#dependencies)
- [Install](#install)
  - [No installation needed - precompiled portable binary](#no-installation-needed---precompiled-portable-binary)
  - [pip](#pip)
- [Running on a sample example](#running-on-a-sample-example)
- [Usage](#usage)

## Dependencies

`make_prg` has two commands: `from_msa` and `update`. The `update` command requires `MAFFT` to be in your `PATH`. It can be installed:
  1. from source: https://mafft.cbrc.jp/alignment/software/;
  2. using `conda`: `conda install -c bioconda mafft`;

## Install

### No installation needed - precompiled portable binary

You can use `make_prg` with no installation at all by simply downloading the precompiled binary, and running it.
In this binary, all libraries are linked statically.

* **Requirements**:
  * `GLIBC >= 2.17` (present on `Ubuntu >= 13.04`, `Debian >= 8.0`, `CentOS >= 7`, `RHEL >= 7.9`,
  `Fedora >= 19`, etc);

* **Download**:
  ```
  wget https://github.com/leoisl/make_prg/releases/download/v0.3.0/make_prg_0.3.0
  ```
* **Running**:
```
chmod +x make_prg_0.3.0
./make_prg_0.3.0 -h
```

* **Credits**:
  * Compilation is done using [PyInstaller](https://github.com/pyinstaller/pyinstaller).

* **Notes**:
  * We provide precompiled binaries for Linux OS only;


### pip

* **Requirements**: `python>=3.7`

* **Installing**:
```sh
pip install git+https://github.com/leoisl/make_prg
```

* **Running**:
```
make_prg -h
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
  
    from_msa     Make PRG from multiple sequence alignment dir
    update       Update PRGs given new sequences output by pandora.

```

#### `from_msa`

```
$ make_prg from_msa --help
usage: make_prg from_msa

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input dir: all files in this will try to be read as the supported alignment_format. If not aligned in fasta alignment_format, use -f to input the alignment_format type
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Output prefix: prefix for the output files
  -t THREADS, --threads THREADS
                        Number of threads
  -f ALIGNMENT_FORMAT, --alignment_format ALIGNMENT_FORMAT
                        Alignment format of MSA, must be a biopython AlignIO input alignment_format. See http://biopython.org/wiki/AlignIO. Default: fasta
  --max_nesting MAX_NESTING
                        Maximum number of levels to use for nesting. Default: 5
  -L MIN_MATCH_LENGTH, --min_match_length MIN_MATCH_LENGTH
                        Minimum number of consecutive characters which must be identical for a match. Default: 7
  -v, --verbose         Increase output verbosity
```

#### `update`

```
$ make_prg update --help
usage: make_prg update

optional arguments:
  -h, --help            show this help message and exit
  -u UPDATE_DS, --update_DS UPDATE_DS
                        Filepath to the update data structures. Should point to a file *.update_DS.
  -d DENOVO_PATHS, --denovo_paths DENOVO_PATHS
                        Filepath containing denovo sequences output by pandora. Should point to a denovo_paths.txt file.
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Output prefix: prefix for the output files
  -t THREADS, --threads THREADS
                        Number of threads
  --mafft MAFFT         Path to MAFFT executable. By default, it is assumed to be on PATH
  --keep_temp           Keep temp files.
  -v, --verbose         Increase output verbosity
```

[pandora]: https://github.com/rmcolq/pandora

