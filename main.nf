include { logSuccess
    logWarning
    logError
    logInfo
    mkdirs
    coloredTitle
    getAbsolute
} from './utils.nf'

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
        Usage: nextflow run papaemmelab/nf-ffperase --step <preprocess|classify|full> [options]

        FFPErase is a tool for pre-processing and classifying FFPE artifact mutations.
        
        For more info: https://github.com/papaemmlab/nf-ffperase
        
        Steps (--step):
            full                [Default] Runs both preprocess + classify.
            preprocess          Compute metrics for features from input vcf.
            classify            Classify variants as real or FFPE artifacts.

        Preprocess Options:
            --vcf               SNV/indel mutations VCF file [required].
            --bam               Input FFPE bam file [required].
            --reference         Reference fasta used to align the FFPE bam [required].
            --outdir            Output location for results [required].
            --coverage          Calculated median coverage [required].
            --medianInsert      Calculated median insert size [required].
            --mutationType      Mutation type, valid choices: "snvs", "indels". [Default: "snvs"]
            --bed               Bedfile path for the regions covered by the bam.
                                [default: assets/gr37.no_mt_unmapped.bed.gz]
            --minBaseq          Minimum BaseQ to assess reads with pileup. [0-60] [default: 20]
            --minDepth          Minimum read depth to assess reads with pileup. [default: 20]
            --minMapq           Minimum MAPQ to assess reads with pileup. [0-60] [default: 20]
            --splitPileup       Number of variants per file for pileup jobs. [default: 50]
            --splitReads        Number of reads to split into picard jobs. [default: 7,500,000]
            --picardMetrics     Output to pre-computed Picard's CollectSequencingArtifactMetrics.

        Classify Options:
            --features          Tsv with preprocessed features. [default: <outdir>/preprocess/features.tsv]
            --model             Path to trained model [required].
            --modelName         Name of the trained model [required].
            --outdir            Output location for results [required].
            --tsv               Tsv that will be used to add annotated columns to the classified output.

    """.stripIndent()
    exit 0
}

def showInfo() {
    logMessage = """\
        ================================================================
        ${coloredTitle('  ')}
        ================================================================

        Documentation @ https://github.com/papaemmlab/nf-ffperase
        Log issues    @ https://github.com/papaemmlab/nf-ffperase/issues

        ----------------------------------------------------------------
        Workflow:
        ----------------------------------------------------------------
        projectDir    : ${workflow.projectDir}
        workDir       : ${workflow.workDir}

        Cmd line:
        \$ ${workflow.commandLine}

        ----------------------------------------------------------------
        Running ${params.step} workflow with run parameters:
        ----------------------------------------------------------------
        step          : ${params.step}
        outdir        : ${params.outdir}
    """

    logMessage += ["preprocess", "full"].contains(params.step) ? (
    """
        ** To Preprocess: **

        bam           : ${params.bam}
        reference     : ${params.reference}
        vcf           : ${params.vcf}
        coverage      : ${params.coverage}
        medianInsert  : ${params.medianInsert}
        mutationType  : ${params.mutationType}
        bed           : ${params.bed}
        picard        : ${params.picard}
        picardMetrics : ${params.picardMetrics ? params.picardMetrics : "''"}
        minMapq       : ${params.minMapq}
        minBaseq      : ${params.minBaseq}
        minDepth      : ${params.minDepth}
        splitReads    : ${params.splitReads}
        splitPileup   : ${params.splitPileup}
    """) : ""
    
    logMessage += ["classify", "full"].contains(params.step) ? (
    """
        ** To Classify: **

        features      : ${params.features ? params.features : "''"}
        model         : ${params.model}
        modelName     : ${params.modelName}
        tsv           : ${new File(params.tsv).name != 'NO_FILE' ? params.tsv : "''"}\
    """) : ""

    logMessage += """
        ----------------------------------------------------------------
    """

    log.info(logMessage.stripIndent())
}

def validateInputs() {

    requiredParams = [
        outdir: true,
    ]
    if (["preprocess", "full"].contains(params.step)) {
        requiredParams.put("bam", true)
        requiredParams.put("reference", true)
        requiredParams.put("vcf", true)
        requiredParams.put("bed", false)
        requiredParams.put("picard", false)
        requiredParams.put("picardMetrics", false)
        requiredParams.put("model", true)
    }
    if (["classify", "full"].contains(params.step)) {
        requiredParams.put("features", false)
        requiredParams.put("model", true)
        requiredParams.put("tsv", false)
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
    if (!validSteps.contains(params.step)) {
        logError """\
            Error: Invalid step '${params.step}'
            Available steps: ${validSteps.join(', ')}
        """.stripIndent()
        exit 1
    }
}

workflow preprocessWorkflow {
    main:
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

        picardInputs = splitBed
            .splitText()
            .map{ line -> line.trim() }
            .combine(inputs.bam)
            .combine(inputs.bai)
            .combine(inputs.reference)
            .combine(inputs.picard)
            .map { nested -> nested.flatten() }

        splitMetrics = picard(picardInputs) | collect

        picardOutput = merge_picard(
            splitMetrics
        )
    }

    // 3. Annotate with Pileup and Picard results
    featuresTsv = annotate_variants(
        pileupOutput,
        picardOutput,
        inputs.reference,
    )

    emit:
    featuresTsv
}

workflow classifyWorkflow {
    take:
    featuresTsv
    
    main:
    inputs = validateInputs()

    classifiedTsv = classify_random_forest(
        featuresTsv,
        inputs.model,
        params.modelName,
        params.mutationType,
        inputs.tsv
    )

    emit:
    classifiedTsv
}

workflow {
    if (params.help) { showHelp() }
    if (params.version) { showVersion() }

    validateSteps()
    showInfo()
    
    def featuresTsv

    if (params.step == "preprocess") {
        featuresTsv = preprocessWorkflow()
    }

    if (params.step == "classify") {
        featuresTsv = inputs.features
        classifyWorkflow(featuresTsv)
    }

    if (params.step == "full") {
        featuresTsv = preprocessWorkflow()
        classifyWorkflow(featuresTsv)
    }
}

workflow.onComplete {
    workflow.success
        ? logSuccess("\nDone! ${coloredTitle()}\u001B[32m ran successfully. See results: ${getAbsolute(params.outdir)}")
        : logError("\nOops .. something went wrong")
}
