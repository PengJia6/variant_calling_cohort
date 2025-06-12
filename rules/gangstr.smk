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
        unpack(get_samples_bam),
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
        bed=lambda wildcards: config["refs"][wildcards.ref_name]["tandem_repeat"]["GangSTR"]
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/gangstr/chroms/{cohort}.{ref_name}.{tech}.gangstr.{chrom}.TR.raw.vcf.gz",
    log:
        config["dir_data"] + "variants_raw/{cohort}/gangstr/logs/{cohort}.{ref_name}.{tech}.{chrom}.gangstr.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/gangstr/logs/{cohort}.{ref_name}.{tech}.{chrom}.gangstr.rtime.tsv"
    threads: get_run_threads("gangstr_call")
    run:
        pref = f"{output}".rstrip(".vcf.gz")
        gangstr = config["software"]["gangstr"]
        bcftools = config["software"]["bcftools"]
        bam_str = ",".join([i for i in input.bams])
        shell(
            "{gangstr} --bam {bam_str} --ref {input.ref} --regions  {input.bed} --out {pref} "
            "--max-proc-read 50000 "
            " 2>{log} 1>{log} "
        )
        shell("{bcftools} view -o {output} -Oz {pref}.vcf")
