include { logSuccess; logWarning; logError; logInfo; mkdirs } from './utils.nf'

include {
    split_pileup
    pileup
    merge_pileup
} from './modules/pileup.nf'

include {
    split_intervals
    picard
    merge_picard
    copy_picard
} from './modules/picard.nf'

include {
    annotate_variants
} from './modules/annotate.nf'

include {
    classify_random_forest
} from './modules/classify.nf'


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
            full              Run 'preprocess' + 'classify'.

        Preprocess Options:
            --vcf               SNV/indel mutations VCF file [required].
            --bam               Input FFPE bam file [required].
            --reference         Reference fasta used to align the FFPE bam [required].
            --outdir            Output location for results [required].
            --bed               Bedfile path for the regions covered by the bam.
                                [default: assets/gr37.no_mt_unmapped.bed.gz]
            --minBaseq          Minimum BaseQ to assess reads with pileup. [0-60] [default: 20]
            --minDepth          Minimum read depth to assess reads with pileup. [default: 20]
            --minMapq           Minimum MAPQ to assess reads with pileup. [0-60] [default: 20]
            --splitReads        Number of reads to split into picard jobs. [default: 7,500,000]
            --coverage          Calculated median coverage [required].
            --medianInsert      Calculated median insert size [required].
            --picardMetrics     Output to pre-computed Picard's CollectSequencingArtifactMetrics.
            --mutationType      Mutation type, valid choices: "snvs", "indels". [Default: "snvs"]

        Classify Options:
            --model             Path to trained model [required].
            --modelName         Name of the trained model [required].
            --outdir            Output location for results [required].

    """.stripIndent()
    exit 0
}

def showInfo() {
    stepInputs = params.step == "preprocess" ? (
     """vcf              : ${params.vcf}
        bed              : ${params.bed}"""
    ) : (
     """model            : ${params.model}"""
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

def validateInputs() {

    requiredParams = [
        bam: true,
        reference: true,
        outdir: true,
    ]
    if (params.step == "preprocess") {
        requiredParams.put("vcf", true)
        requiredParams.put("bed", false)
        requiredParams.put("picard", false)
        requiredParams.put("picardMetrics", false)
        requiredParams.put("model", true)
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
    def validSteps = ['preprocess', 'classify', 'full']
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

workflow {
    if (params.help) { showHelp() }
    if (params.version) { showVersion() }

    validateSteps()
    showInfo()
    
    switch (params.step) {
        case 'preprocess':
            preprocessWorkflow()
            break

        case 'classify':
            classifyWorkflow()
            break

        case 'full':
            preprocessWorkflow()
            classifyWorkflow()
            break
    }
}

workflow preprocessWorkflow {
    // Check Inputs
    if (!params.coverage) {
        logError "Error: --coverage is required."
        exit 1
    }
    if (!params.medianInsert) {
        logError "Error: --medianInsert is required."
        exit 1
    }
    def validMutationTypes = ["snvs", "indels"]
    if (!validMutationTypes.contains(params.mutationType)) {
        logError """\
            Error: Invalid Mutation Type: '${params.mutationType}.'
            Valid choices are: ${validMutationTypes.join(', ')}.
        """.stripIndent()
        exit 1
    }

    inputs = validateInputs()

    // 1. Pileup Mutations
    splitVcfs = split_pileup(
        inputs.vcf
    ) | flatten

    pileupInputs = splitVcfs
        .combine(inputs.bam)
        .combine(inputs.bai)
        .combine(inputs.reference)
        .map { nested -> nested.flatten() }
    pileupVcfs = pileup(pileupInputs) | collect

    pileupOutput = merge_pileup(pileupVcfs)

    // 2. Get Metrics from Picard
    if (inputs.picardMetrics) {
        // Read from pre-computed metrics
        picardOutput = copy_picard(
            inputs.picardMetrics
        )
    } else {
        // Compute new metrics
        splitBed = split_intervals(
            inputs.bam,
            inputs.bai,
            inputs.bed,
        )
        splitMetrics = picard(
            splitBed.splitText().map { line -> line.trim() },
            inputs.bam,
            inputs.bai,
            inputs.reference,
            inputs.picard,
        )
        picardOutput = merge_picard(
            splitMetrics
        )
    }

    // 3. Annotate with Pileup and Picard results
    annotatedVariants = annotate_variants(
        pileupOutput,
        picardOutput,
        inputs.reference,
        params.mutationType == "indels"
    )

    // 4. Classify
    classify_random_forest(
        annotatedVariants,
        inputs.model,
        params.modelName,
        params.mutationType
    )
}

workflow classifyWorkflow {
    inputs = validateParams()

    // classify(
    //     inputs.tsv,
    //     inputs.model,
    //     params.model_name,
    //     params.mutationType,
    // )
}

workflow.onComplete {
    workflow.success 
        ? logSuccess("\nDone! FFPErase ran successfully. See the results in: ${params.outdir}") 
        : logError("\nOops .. something went wrong")
}