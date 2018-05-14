# make_prg
Code to create a PRG for input to Pandora from a Multiple Sequence Alignment file.

__Requirements__
Expects you to have python3 and nextflow installed and in your path.

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


