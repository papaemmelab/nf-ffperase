include { logSuccess; logWarning; logError; logInfo; mkdirs } from './utils.nf'

include { split_pileup } from './modules/pileup.nf'
include { pileup } from './modules/pileup.nf'
include { merge_pileup } from './modules/pileup.nf'

// params.outdir = "${workflow.projectDir}/results"

def showVersion() {
    version = "v0.1.0"
    
    log.info """\
        Repo: https://github.com/papaemmlab/nf-ffperase
        FFPErase version: ${version}
    """.stripIndent()
    exit 0
}

def showHelp() {
    logInfo """\
        Usage: nextflow run papaemmelab/nf-ffperase --step <preprocess|classify> [options]

        FFPErase is a tool for pre-processing and classifying FFPE artifact mutations.
        
        For more info: https://github.com/papaemmlab/nf-ffperase

        Steps:
            preprocess        Compute metrics from input vcf.
            classify          Classify variants as real or FFPE artifacts.

        Preprocess Options:
            --vcf               SNV/indel mutations VCF file [required].
            --bam               Input FFPE bam file [required].
            --reference         Reference fasta used to align --ffpe-bam [required].
            --outdir            Output location for results [required].

        Classify Options:
            --model             Path to trained model [required].
            --bam               Input FFPE bam file [required].
            --reference         Reference fasta used to align --ffpe-bam [required].
            --outdir            Output location for results [required].

    """.stripIndent()
    exit 0
}

def showInfo() {
    stepInputs = params.step == "preprocess" ? (
        "vcf              : ${params.vcf}"
    ): (
        "model            : ${params.model}"
    )

    log.info """\
        ===================================================================
        F F P E R A S E
        ===================================================================

        Documentation   @ https://github.com/papaemmlab/nf-ffperase
        Log issues      @ https://github.com/papaemmlab/nf-ffperase/issues

        -------------------------------------------------------------------
        Workflow run parameters
        -------------------------------------------------------------------
        step             : ${params.step}
        workDir          : ${workflow.workDir}
        -------------------------------------------------------------------
        ${stepInputs}
        bam              : ${params.bam}
        reference        : ${params.reference}
        outdir           : ${params.outdir}
        -------------------------------------------------------------------
        
        Running ${params.step} workflow...
    """.stripIndent()
}

def validateParams() {

    requiredParams = [
        bam: true,
        reference: true,
        outdir: true,
    ]
    if (params.step == "preprocess") {
        requiredParams.put("vcf", true)
    } else {
        requiredParams.put("model", true)
    }

    def channels = [:]
    requiredParams.each { key, isRequired ->
        // Validate param is provided
        def paramValue = params."${key}"
        if (!paramValue && isRequired) {
            logError "Error: --${key} is required."
        }

        // If provided validate path exists
        if (paramValue) {
            channels[key] = channel.fromPath(paramValue)
                .ifEmpty("Error: No valid file found for ${key} at ${paramValue}")
        } else {
            channels[key] = null
        }

        // check for bam index
        if (key == "bam") {
            channels["bai"] = channel.fromPath("${paramValue}.bai")
                .ifEmpty("Error: No valid file found for bai at ${paramValue}.bai")
        }
    }
    return channels
}

def validateSteps() {
    def validSteps = ['preprocess', 'classify']
    if (!params.step) {
        showHelp()
    } else if (!validSteps.contains(params.step)) {
        logError """\
            Error: Invalid step '${params.step}'
            Available steps: ${validSteps.join(', ')}
        """.stripIndent()
        exit 1
    }
}

def createOutdirs() {
    if (!params.outdir) {
        params.outdir = "${workflow.projectDir}/results"
    }
    
    outdirPreprocess = "${params.outdir}/preprocessed_results"
    outdirClassify = "${params.outdir}/classification_results"
    
    switch (params.step) {
        case 'preprocess':
            mkdirs(outdirPreprocess)
            break

        case 'classify':
            mkdirs(outdirClassify)
            break
        
        case 'both':
            mkdirs(outdirPreprocess)
            mkdirs(outdirClassify)
            break
    }
}

workflow {
    if (params.help) { showHelp() }
    if (params.version) { showVersion() }

    validateSteps()
    createOutdirs()
    showInfo()
    
    switch (params.step) {
        case 'preprocess':
            preprocessWorkflow()
            break

        case 'classify':
            classifyWorkflow()
            break
    }
}

workflow preprocessWorkflow {
    inputs = validateParams()

    // Pileup Mutations
    splitVcfs = split_pileup(
        inputs.vcf, 
        inputs.outdir
    ) | flatten

    pileupInputs = splitVcfs
        .combine(inputs.bam)
        .combine(inputs.bai)
        .combine(inputs.reference)
        .combine(inputs.outdir)
        .map { nested -> nested.flatten() }
    pileupVcfs = pileup(pileupInputs) | collect

    mergedPileup = merge_pileup(pileupVcfs, inputs.outdir)

    
}

workflow classifyWorkflow {
    inputs = validateParams()

    classify(
        params.bam,
        params.reference,
        params.model,
        params.outdir,
        params.model_name
    )
}

process classify {
    input:
    path bam
    path reference
    path model
    path outdir
    // val modelName

    output:
    path "${outdir}/classification_results"

    script:
    """
    echo "Classifying FFPE BAM file: $bam"
    echo "Using reference genome: $reference"
    echo "Using model: $model ($modelName)"
    echo "Output will be stored in ${outdir}/classification_results"
    """
}

workflow.onComplete {
    workflow.success 
        ? logSuccess("\nDone! FFPErase ran successfully. See the results in: ${params.outdir}") 
        : logError("\nOops .. something went wrong")
}