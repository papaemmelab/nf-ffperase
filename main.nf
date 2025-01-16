// Import modules
include { logSuccess; logWarning; logError; logInfo } from './utils.nf'
// include {  split_pileup, pileup, merge_pileup } from './modules/pileup.nf'

params.outdir = "${workflow.projectDir}/results"

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
            --outdir            Output location for results [default: ./results].

        Classify Options:
            --model             Path to trained model [required].
            --bam               Input FFPE bam file [required].
            --reference         Reference fasta used to align --ffpe-bam [required].
            --outdir            Output location for results [default: ./results].
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
        version          : ${params.version}
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
        outdir: false
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
        log.info "Checking ${key}: ${paramValue}"
        
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
    }
}

workflow preprocessWorkflow {
    inputs = validateParams()

    preprocess(
        inputs.bam,
        inputs.reference,
        inputs.vcf,
        params.outdir
    )
}

workflow classifyWorkflow {
    classify(
        params.bam,
        params.reference,
        params.model,
        params.outdir,
        params.model_name
    )
}
    
process preprocess {
    input:
    path vcf
    path bam
    path reference
    path outdir

    output:
    path "${outdir}/preprocessed_results"

    script:
    """
    echo "Preprocessing FFPE BAM file: $bam"
    echo "Using reference genome: $reference"
    echo "Input VCF: $vcf"
    echo "Output will be stored in ${outdir}/preprocessed_results"

    mkdir -p ${outdir}/preprocessed_results
    """
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


// Define the preprocess workflow
// def runPreprocess() {
//     println "Running preprocessing step..."
//     process preprocess {
//         input:
//         path ffpeBam from params.bam
//         path reference from params.reference
//         path vcf from params.vcf

//         output:
//         path "${params.outdir}/preprocessed_results"

//         script:
//         """
//         # Example preprocessing script
//         echo "Preprocessing FFPE BAM file: $ffpeBam"
//         echo "Using reference genome: $reference"
//         echo "Input VCF: $vcf"
//         echo "Output will be stored in ${params.outdir}/preprocessed_results"
//         """
//     }
// }

// // Define the classify workflow
// def runClassify() {
//     println "Running classification step..."
//     process classify {
//         input:
//         path ffpeBam from params.bam
//         path reference from params.reference
//         path modelPath from params.model
//         val modelName from params.model_name

//         output:
//         path "${params.outdir}/classification_results"

//         script:
//         """
//         # Example classification script
//         echo "Classifying FFPE BAM file: $ffpeBam"
//         echo "Using reference genome: $reference"
//         echo "Using model: $modelPath ($modelName)"
//         echo "Output will be stored in ${params.outdir}/classification_results"
//         """
//     }
// }