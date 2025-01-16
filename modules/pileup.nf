process split_pileup {
    input:
    path vcf

    output:
    path "splits/split_*.vcf", emit: split_vcf_files

    script:
    """
    #!/usr/bin/env python
    ANNOT_LIMIT = 1  # Adjust this limit as needed

    vcf = '${vcf}'
    split_dir = 'splits'
    import os
    os.makedirs(split_dir, exist_ok=True)

    with open(vcf, 'r', encoding='utf-8') as f_in:
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
    input:
    tuple path(split_vcf), path(input_bam), path(input_bai), path(ref_genome)

    output:
    path "pileup/pileup_*.txt", emit: pileup_files

    script:
    """
    mkdir -p pileup
    BASENAME_VCF=\$(basename ${split_vcf} .vcf)
    PILEUP_FILE="pileup/pileup_\${BASENAME_VCF}.txt"
    annotate_w_pileup ${input_bam} ${ref_genome} ${split_vcf} \$PILEUP_FILE --mapq 0 --baseq 0 --depth 0
    """
}

process merge_pileup {
    input:
    path pileup_files

    output:
    path "merged_pileup.txt", emit: merged_pileup

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
