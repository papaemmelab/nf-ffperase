// Import modules
include { split_pileup } from './modules/pileup.nf'
include { pileup } from './modules/pileup.nf'
include { merge_pileup } from './modules/pileup.nf'


// Header
log.info """\

===================================================================
F F P E R A S E
===================================================================

Documentation	@ https://github.com/papaemmlab/nf-ffperase

Log issues		@ https://github.com/papaemmlab/nf-ffperase/issues

===================================================================
Workflow run parameters
===================================================================

 version          : ${params.version}
 inputVcf         : ${params.inputVcf}
 inputBam         : ${params.inputBam}
 reference        : ${params.reference}
 outDir           : ${params.outDir}
 workDir          : ${workflow.workDir}

===================================================================
 """

// Help function
def helpMessage() {
log.info"""

Usage:  nextflow run main.nf --input samplesheet.tsv --ref reference.fasta

Required Arguments:

	--input			    Full path and name of sample input file (tsv format).

	--ref			    Full path and name of reference genome (fasta format).

	--outDir		    Full path and name of results directory.

Optional Arguments:

	--version 		    Show version and exit.

""".stripIndent()
}

workflow {

    if (params.help == true) {
        helpMessage()
        exit 1
    } else {
        // Create channels for the inputs
        input_vcf = Channel.fromPath(params.inputVcf)
        input_bam = file(params.inputBam)
        input_bai = file(params.inputBam + ".bai")
        ref_genome = file(params.reference)

        // Split VCF file into multiple files
        split_vcf_files = split_pileup(input_vcf)

		// Package each split VCF file into a tuple with the BAM, BAI, and reference genome
        split_vcf_files
            .flatten()
            .map { split_vcf -> tuple(split_vcf, input_bam, input_bai, ref_genome) }
            .set { pileup_input }

        // pileup_input.view()
		pileup(pileup_input)

        // Merge the pileup output files
        merge_pileup(pileup_output.collect())
    }
}

workflow.onComplete {

summary = """
===================================================================
Workflow execution summary
===================================================================

Duration		: ${workflow.duration}
Success			: ${workflow.success}
workDir			: ${workflow.workDir}
Exit status	: ${workflow.exitStatus}
outDir			: ${params.outDir}

===================================================================

"""

println summary

}