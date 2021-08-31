# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `-S`, `--seqid` option to name the PRG sequence, which by default uses the file name.
- `-N` shortcut for max nesting
- `-L` shortcut for min match length
- `--log` to enable specifying log file should go to path. Default behaviour is now that
  log goes to stderr by default
- `--summary` option to specify that user wants a summary file, rather than the default
  behaviour where they get no choice
- `-O`, `--output-type` option to specify what output files are required. Defaults to
  all

### Changed
- `--prefix` CLI parameter of `from_msa` subcommand removed in favor of CLI parameters `--outdir`
   and `--prg_name`, with sensible defaults (current working directory and MSA file name stem respectively).
   This allows finer control over where to place output files.
- Output files no longer contain 'max_nesting' and 'min_match_length' in their names; these appear in the log files,
  and in the `.prg` fasta header.
- `.bin` file now stores even integer markers at site ends; this is the format used by gramtools.
- Summary file not written by default

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


[Unreleased]: https://github.com/rmcolq/make_prg/compare/v0.1.1...HEAD

[0.1.1]: https://github.com/rmcolq/make_prg/releases/v0.1.1
[0.1.0]: https://github.com/rmcolq/make_prg/releases/v0.1.0
