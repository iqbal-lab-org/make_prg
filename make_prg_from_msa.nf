#!/usr/bin/env nextflow

params.help = false
params.ref_dir = ""
params.output_dir = ""
params.max_forks = 100

if (params.help){
    log.info"""
	Pipeline to create a Pandora-style Pangenome PRG file from a directory of MSA
	e.g. downloaded for a species of interest from the panX website

	Required arguments:
	  --ref_dir DIRECTORY Directory containing input MSA in fasta format
	  --output_dir DIRECTORY Directory where the final pangenome PRG fasta should be output
	Optional arguments:
	  max_forks INT Limit number of MSA that can be processed concurrently

    """.stripIndent()
    exit 0
}

if (!ref_dir.isDirectory()) {
    exit 1, "Directory containing input MSA in fasta format not found: ${params.ref_dir} -- aborting"
}

if (!output_dir.isDirectory()) {
    exit 1, "Directory for output not found: ${params.output_dir} -- aborting"
}
