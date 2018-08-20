# make_prg
Code to create a PRG for input to Pandora (https://github.com/rmcolq/pandora) from a Multiple Sequence Alignment file.

__Requirements__
Expects you to have python3 and nextflow installed and in your path and a config file for nextflow set up if you are working on a cluster. Can run on python2.7+ if the command in nextflow file is edited.
Note that nextflow does not play nicely when files are in mounted or shared folders.

__Usage__

    Usage: nextflow run make_prg_nexflow.nf <arguments>
  
    Required arguments:
      --tsv_in  FILENAME  An index file of MSA to build PRGs of
      --pipeline_root DIRECTORY Absolute path to make_prg
    
    Optional arguments:
  
__Download__
```
git clone https://github.com/rmcolq/make_prg.git
cd make_prg
pip3 install -r requirements.txt
pytest 
```

__Input__
Multiple Sequence Alignment files for genes/dna sequences for which we want PRGs, and an tab-separated index of these in the form:
```
sample_id       infile
GC0000001   /absolute/path/to/GC0000001_na_aln.fa.gz
GC0000002   /absolute/path/to/GC0000002_na_aln.fa
```

__Changing parameters__

There are some parameters at the top of the nextflow file which could be changed but which I have not made command line parameters:
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
