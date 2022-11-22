# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2021-26-11

### Added
- `make_prg update` command, that updates PRGs without requiring to rebuild MSAs and the PRG itself from scratch;
- Trace (`-vv`) logging level, to track make_prg behaviour (intended for developers only);
- Multithreading support (`-t` parameter);
- `--output_graphs` parameter, that outputs the recursive tree graphical representation (intended for developers only);
- A sample example;
- 260 new tests (from 116 to 376 total tests), with test coverage of around 99% in non-argument parsing code;

### Changed
- `make_prg from_msa`: input can now be a single file or a directory. If it is a single file,
a `<prefix>.prg.bin`, a `<prefix>.prg.fa`, a `<prefix>.prg.gfa` and a `<prefix>.update_DS.zip`
files are created. If it is a directory, all files in the directory are scanned and the same
execution for a single file is done for each input file found. The output files are a collection of the single-input
execution: a `<prefix>.prg.bin.zip` file will contain a collection of `.prg.bin` files, similar to
`<prefix>.prg.gfa.zip` and `<prefix>.update_DS.zip`; `<prefix>.prg.fa` will be a multi fasta;
- The recursive clustering and collapse algorithm is now explicitly represented as a tree with internal data
structures that remember the multiple sequence subalignment at any point of the recursion, as well as several other
internal data, allowing the serialization and deserialization of the recursion tree at any point. Thus updates can
be done avoiding any recomputation by firstly saving the state of the recursion tree to disk, and then loading this
recursion tree, adding denovo sequences to some specific nodes, and triggering recomputation of the modified nodes.
Any preorder traversal of the recursion tree yields the same order of recursive calls of the previous algorithm,
thus allowing us to translate the algorithms in the previous version as preorder traversals with custom visit
operations.

### Fixed
- Several minor bugs;
- Heavy refactoring of almost the whole codebase;

## [0.2.0] - 2021-10-27

### Added
- New CLI:
    - `-S`, `--seqid` option to name the PRG sequence, which by default uses the file name.
    - `-N` shortcut for max nesting
    - `-L` shortcut for min match length
    - `--log` to enable specifying log file should go to path. Default behaviour is now that
      log goes to stderr by default
    - `-O`, `--output-type` option to specify what output files are required. Defaults to
      all
- Integration test for nested variation

### Removed
- summary file


### Changed
- New CLI:
    - `--prefix` CLI parameter of `from_msa` subcommand removed in favor of CLI parameters `--outdir`
       and `--prg_name`, with sensible defaults (current working directory and MSA file name stem respectively).
       This allows finer control over where to place output files.
- Output files:
    - No longer contain 'max_nesting' and 'min_match_length' in their names; these appear in the log files,
      and in the `.prg` fasta header.
    - `.bin` file now stores even integer markers at site ends; this is the format used by gramtools.
    - Summary file not written by default

### Fixed
- Bugfix ([#27][27])- Part 1: clustering termination was not properly detected, causing spurious 
  site production
- Bugfix ([#27][27])- Part 2: added criteria for not clustering, to reduce construction of 
  ambiguous graphs
  

## [0.1.1] - 2021-01-22
### Added
- Dockerfile
- `-V` option to get version

### Changed
- A test that was clustering all unique 5-mers was reduced to all 4-mers as the memory
  usage of all 5-mers was causing a segfault when trying to run the tests during the
  docker image build.

### Removed
- Singularity file as it is redundant with the new Dockerfile (that will be hosted on
  quay.io)
- `scipy` dependency. We never actually explicitly use `scipy`.

## [0.1.0] - 2021-01-20
### Added
- This CHANGELOG file to hopefully serve as an evolving example of a standardized open
  source project CHANGELOG.


[Unreleased]: https://github.com/iqbal-lab-org/make_prg/compare/v0.1.0...HEAD

[1.0.0]: https://github.com/iqbal-lab-org/make_prg/releases/tag/0.1.0
[0.2.0]: https://github.com/iqbal-lab-org/make_prg/releases/tag/0.2.0
[0.1.1]: https://github.com/iqbal-lab-org/make_prg/releases/tag/0.1.1
[0.1.0]: https://github.com/iqbal-lab-org/make_prg/releases/tag/0.1.0

[27]: https://github.com/iqbal-lab-org/make_prg/issues/27
