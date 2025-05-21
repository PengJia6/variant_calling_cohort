rule manta_conf:
    input:
        unpack(get_samples_bam),
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        prefix=config["dir_data"] + "variants_raw/{cohort}/manta/chroms/{cohort}.{ref_name}.{tech}.{chrom}.manta/runWorkflow.py"
    params:
        extra="",
    # dp=5
    log:
        config["dir_data"] + "variants_raw/{cohort}/manta/logs/{cohort}.{ref_name}.{tech}.{chrom}.manta_conf.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/manta/logs/{cohort}.{ref_name}.{tech}.{chrom}.manta_conf.rtime.tsv"
    threads: get_run_threads("manta_call")
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
        prefix = str(output.prefix)[:-15]
        shell("mkdri -p {prefix}")
        manta = config["software"]["manta"]
        bgzip = config["software"]["bgzip"]
        tabix = config["software"]["tabix"]

        fasta = pysam.FastaFile(f"{input.ref}")
        chrom_list = {j: i for i, j in zip(fasta.references,fasta.lengths)}
        fasta.close()
        bed_chrom = open(f"{prefix}/{wildcards.chrom}.bed","w")
        bed_chrom.write(f"{wildcards.chrom}\t0\t{chrom_list[wildcards.chrom]}\n")
        bed_chrom.close()
        shell("{bzip} {prefix}/{wildcards.chrom}.bed")
        shell("{tabix} {prefix}/{wildcards.chrom}.bed.gz")
        bams_str=" ".join([f"--bam {i}" for i in input.bam])
        shell("{manta} --bam {input.bams} --referenceFasta {input.ref} "
              "--runDir {prefix} --callRegions {prefix}/{wildcards.chrom}.bed.gz 2>{log} 1>{log}")

rule manta:
    input:
        prefix=rules.manta_conf.output.prefix
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/manta/chroms/{cohort}.{ref_name}.{tech}.{chrom}.manta.raw.vcf.gz"
    log:
        config["dir_data"] + "variants_raw/{cohort}/manta/logs/{cohort}.{ref_name}.{tech}.{chrom}.manta.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/manta/logs/{cohort}.{ref_name}.{tech}.{chrom}.manta.rtime.tsv"

    threads: get_run_threads("manta")
    run:
        prefix = str(input.prefix)[:-15]
        tmp_vcf = f"{prefix}/results/variants/candidateSV.vcf.gz",
        configmanta = config["software"]["manta"]
        manta_pre = "/".join(f"{configmanta}".split("/")[:-1])
        shell("export PATH={manta_pre}:$PATH && "
              "chmod +x {input.prefix} && "
              "{input.prefix} -j {threads} 2>{log} 1>{log}")
        shell("ln {tmp_vcf} {output.vcfgz}")
        shell("ln {tmp_vcf}.tbi {output.vcfgz}.tbi")
