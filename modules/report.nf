process PLOT_REPORT {
    publishDir "${params.outdir}/plots", mode: "copy"

    container "papaemmelab/pycirclize:1.9.1"
    
    input:
    path classifiedTsv
    val mutationType
    
    output:
    path "distributions.png", emit: distributionPlot
    path "circos.png", emit: circosPlot
    
    script:
    template 'plot_report.py'
}
