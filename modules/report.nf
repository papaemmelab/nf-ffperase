process PLOT_REPORT {
    publishDir "${params.outdir}/plots", mode: "copy"

    container "papaemmelab/pycirclize:1.9.1"
    
    input:
    path classifiedTsv
    val mutationType
    path chrom_bed from file("${workflow.projectDir}/assets/hg19.chrom.bed")
    path cytoBand from file("${workflow.projectDir}/assets/cytoBand.txt")
    
    output:
    path "distributions.png", emit: distributionPlot
    path "circos.png", emit: circosPlot
    
    script:
    template 'plot_report.py'
}
