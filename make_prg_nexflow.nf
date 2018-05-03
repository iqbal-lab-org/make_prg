params.tsv_in = ""
params.testing = false
params.pipeline_root = ""
params.max_nesting = 10
params.min_match_length = 7
params.alignment_format = "fasta"
params.final_outdir = "."

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
    errorStrategy {task.attempt < 4 ? 'retry' : 'terminate'}
    memory {params.testing ? '0.5 GB' : 1.GB * task.attempt * task.attempt}
    maxRetries 4

    input:
    val tsv_fields from split_tsv

    output:
    set(val("${tsv_fields['sample_id']}"), file("${tsv_fields['sample_id']}.max_nest${params.max_nesting}.min_match${params.min_match_length}.prg")) into make_prg_out

    """
    python ${params.pipeline_root}/make_prg_from_msa.py ${tsv_fields["infile"]} --prefix ${tsv_fields['sample_id']} --max_nesting ${params.max_nesting} --min_match_length ${params.min_match_length} --alignment_format ${params.alignment_format} -v
    """
}

process make_fasta {
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
