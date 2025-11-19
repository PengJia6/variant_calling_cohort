# ======================================================================================================================
# Project: variant_calling_cohort
# Script : trgt.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: trgt call
# ======================================================================================================================
rule gangstr_call:
    input:
        # unpack(get_samples_bam),
        bam= lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
        bed=lambda wildcards: config["refs"][wildcards.ref_name]["tandem_repeat"]["GangSTR"]
    output:
        # protected(config["dir_data"] + "variants_raw/{cohort}/pindel/chroms/{cohort}.{ref_name}.{tech}.{chrom}.pindel")
        # vcfgz=config["dir_data"] + "variants_raw/{cohort}/hipstr/chroms/{cohort}.{ref_name}.{tech}.hipstr.{chrom}.TR.raw.vcf.gz",
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/gangstr/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.gangstr.TR.raw.vcf.gz",


    log:
        config["dir_data"] + "variants_raw/{cohort}/gangstr/logs/samples/{cohort}.{sample}.{ref_name}.{tech}.{chrom}.gangstr.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/gangstr/logs/samples/{cohort}.{sample}.{ref_name}.{tech}.{chrom}.gangstr.rtime.tsv"
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    # output:
        # vcfgz=config["dir_data"] + "variants_raw/{cohort}/gangstr/chroms/{cohort}.{ref_name}.{tech}.gangstr.{chrom}.TR.raw.vcf.gz",
    # log:
        # config["dir_data"] + "variants_raw/{cohort}/gangstr/logs/{cohort}.{ref_name}.{tech}.{chrom}.gangstr.log"
    # benchmark:
        # config["dir_data"] + "variants_raw/{cohort}/gangstr/logs/{cohort}.{ref_name}.{tech}.{chrom}.gangstr.rtime.tsv"
    threads: get_run_threads("gangstr_call")
    run:
        pref = f"{output}".rstrip(".vcf.gz")
        gangstr = config["software"]["gangstr"]
        bcftools = config["software"]["bcftools"]
        # bam_str = ",".join([i for i in input.bams])
        shell(
            "{gangstr} --bam {input.bam} --ref {input.ref} --regions  {input.bed} --out {pref} "
            "--max-proc-read 3000  --chrom {wildcards.chrom} "
            " 2>{log} 1>{log} "
        )
        shell("{bcftools} view -o {output} -Oz {pref}.vcf")
