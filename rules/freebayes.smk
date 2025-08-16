rule freebayes_call:
    input:
        unpack(get_samples_bam),
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        vcfgz = config["dir_data"] + "variants_raw/{cohort}/freebayes/chroms/{cohort}.{ref_name}.{tech}.freebayes.{chrom}.SNVIndel.raw.vcf.gz"
    params:
        extra="",
    # dp=5
    log:
        config["dir_data"] + "variants_raw/{cohort}/freebayes/logs/{cohort}.{ref_name}.{tech}.{chrom}.freebayes_call.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/freebayes/logs/{cohort}.{ref_name}.{tech}.{chrom}.freebayes_call.rtime.tsv"
    threads: get_run_threads("freebayes_call")
    # input:
    #     bam=config["dir_aligned_reads"] + "{prefix}.{ref_name}.{suffix}.bam",
    #     bai=config["dir_aligned_reads"] + "{prefix}.{ref_name}.{suffix}.bam.bai",
    #     ref=config["dir_ref"] + "{ref_name}.fasta",
    # output:
    #     prefix=config["dir_variants"] + "{prefix}.{ref_name}.{suffix}.manta/runWorkflow.py",
    # log:
    #     config["dir_variants"] + "{prefix}.{ref_name}.{suffix}.manta_conf.log",
    # # benchmark:
    # #          config["dir_logs"] + "dv/{sample}/{sample}.{prefix}.dv.tsv"
    # benchmark:
    #     config["dir_variants"] + "{prefix}.{ref_name}.{suffix}.manta_conf.rtime.tsv",
    #
    # threads: get_run_threads("manta_conf")
    run:
        import pysam
        freebayes = config["software"]["freebayes"]
        bcftools = config["software"]["bcftools"]

        shell("{freebayes}  -r {wildcards.chrom} -f {input.ref} {input.bams}  2>>{log}|"
              "{bcftools} view  --threads {threads} -Oz -o {output.vcfgz} 1>>{log} 2>>{log} ")


