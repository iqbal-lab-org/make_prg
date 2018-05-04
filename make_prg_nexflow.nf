params.tsv_in = ""
params.testing = false
params.pipeline_root = ""
params.max_nesting = 10
params.min_match_length = 7
params.alignment_format = "fasta"
params.final_outdir = "."
params.max_forks_make_prg = 100
params.max_forks_make_fasta = 100

input_tsv = file(params.tsv_in).toAbsolutePath()
if (!input_tsv.isFile()) {
    exit 1, "Directory containing input MSA in fasta format not found: ${params.tsv_in} -- aborting"
}

final_outdir = file(params.final_outdir).toAbsolutePath()
if (!final_outdir.exists()) {
    exit 1, "Output directory not found: ${params.final_outdir} -- aborting"
}

split_tsv = Channel.from(input_tsv).splitCsv(header: true, sep:'\t')

process make_prg {
    maxForks params.max_forks_make_prg
    errorStrategy {task.attempt < 5 ? 'retry' : 'terminate'}
    memory {params.testing ? '0.5 GB' : 0.5.GB * task.attempt * task.attempt}
    maxRetries 5

    input:
    val tsv_fields from split_tsv

    output:
    set(val("${tsv_fields['sample_id']}"), file("${tsv_fields['sample_id']}.max_nest${params.max_nesting}.min_match${params.min_match_length}.prg")) into make_prg_out

    """
    python3 ${params.pipeline_root}/make_prg_from_msa.py ${tsv_fields["infile"]} --max_nesting ${params.max_nesting} --alignment_format ${params.alignment_format} --min_match_length ${params.min_match_length} --prefix ${tsv_fields['sample_id']} -v
    """
}

process make_fasta {
    maxForks params.max_forks_make_fasta
    input: 
    set val(id), file(prg_file) from make_prg_out

    output:
    file "${id}.fa" into results
    
    """
    echo ">${id}" >> "${id}.fa"
    cat "${prg_file}" >> "${id}.fa"
    echo "" >> "${id}.fa"
    """ 
}

results.collectFile(name: final_outdir/'pangenome_PRG.fa')
