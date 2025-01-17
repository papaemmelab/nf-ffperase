process split_pileup {
    publishDir "${outdir}/preprocessed_results", mode: "copy"

    input:
    path vcf
    path outdir

    output:
    path "splits/split_*.vcf", emit: splitVcfs

    script:
    """
    #!/usr/bin/env python
    import os
    
    ANNOT_LIMIT = 1  # Adjust this limit as needed

    split_dir = 'splits'
    os.makedirs(split_dir, exist_ok=True)

    with open('${vcf}', 'r', encoding='utf-8') as f_in:
        lines = list(f_in.readlines())
    header = [l for l in lines if l.startswith('#')]
    records = [l for l in lines if not l.startswith('#')]

    for ix, region in enumerate(range(0, len(records), ANNOT_LIMIT)):
        end = min(len(records), region + ANNOT_LIMIT)
        output_vcf_path = f'{split_dir}/split_{ix}.vcf'
        with open(output_vcf_path, 'w', encoding='utf-8') as f_out:
            f_out.write(''.join(header))
            f_out.write(''.join(records[region:end]))
    """.stripIndent()
}

process pileup {
    publishDir "${outdir}/preprocessed_results", mode: "copy"
    container "papaemmelab/annotate-w-pileup:v0.1.0"
    
    input:
    path splitVcf
    path bam 
    path bai
    path reference

    output:
    path "pileup/pileup_*.txt", emit: pileupVcfs

    script:
    """
    OUTDIR=\$PWD/pileup
    mkdir -p \$OUTDIR
    
    BASENAME_VCF=\$(basename ${splitVcf} .vcf)
    PILEUP_FILE="\$OUTDIR/pileup_\${BASENAME_VCF}.txt"
    
    annotate_w_pileup \
        \$PWD/${bam} \
        \$PWD/${reference} \
        \$PWD/${splitVcf} \
        \$PILEUP_FILE \
        --mapq 0 \
        --baseq 0 \
        --depth 0
    """
}

process merge_pileup {
    publishDir "${outdir}/preprocessed_results", mode: "copy"

    input:
    path pileupVcf

    output:
    path "merged_pileup.txt", emit: mergedPileup

    script:
    """
    # Get the header from the first file
    head -n 1 ${pileup_files[0]} > merged_pileup.txt

    # Concatenate all files, skipping headers from subsequent files
    for f in ${pileup_files[*]}; do
        tail -n +2 "\$f" >> merged_pileup.txt
    done
    """
}
