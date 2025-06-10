# ======================================================================================================================
# Project: variant_calling_cohort
# Script : pindel.smk
# Author : Peng Jia
# Date   :  2024/11/18
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: TODO
# ======================================================================================================================
import random

import numpy as np
import pysam


rule hipstr_call:
    input:
        unpack(get_samples_bam),
        bed=lambda wildcards: config["refs"][wildcards.ref_name]["tandem_repeat"]["HipSTR"]
    # exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
    output:
        # protected(config["dir_data"] + "variants_raw/{cohort}/pindel/chroms/{cohort}.{ref_name}.{tech}.{chrom}.pindel")
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/hipstr/chroms/{cohort}.{ref_name}.{tech}.hipstr.{chrom}.TR.raw.vcf.gz",
    params:
        extra="",
        dp=5
    log:
        config["dir_data"] + "variants_raw/{cohort}/hipstr/logs/{cohort}.{ref_name}.{tech}.{chrom}.hipstr.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/hipstr/logs/{cohort}.{ref_name}.{tech}.{chrom}.hipstr.rtime.tsv"
    threads: get_run_threads("hipstr_call")
    run:
        hipstr = config["software"]["hipstr"]
        bam_str = ",".join([i for i in input.bams])
        shell(
            "{hipstr} --bams {bam_str} --fasta {input.ref} --regions  {input.bed} --str-vcf {output} --chrom {wildcards.chrom} "
            "--output-filters --use-unpaired --no-rmdup --max-str-len 10000 --max-flank-indel 1 >{log} 1>{log} "
        )


# def get_vcfs_for_pindel_merge_vcf(wildcards):
#     gvcfs = []
#     vcfs = []
#     # for cohort, cohort_info in.items():
#     for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
#         # for sample_name, info_t in config["samples_info"][f"{wildcards.cohort}"].items():
#         #     for sample_type, info_tt in info_t["path"].items():
#         # print(wildcards.cohort,sample_type,sample_name,wildcards.this_prefix)
#         # for tech, info in info_tt.items():
#         vcfs.append(config["dir_data"] + f"{wildcards.cohort}/deepvariant/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.deepvariant.SNVIndel.vcf.gz")
#         gvcfs.append(config["dir_data"] + f"{wildcards.cohort}/deepvariant/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.deepvariant.SNVIndel.g.vcf.gz")
#     return {"gvcf_gz": gvcfs, "vcf_gz": vcfs}
# rule pindel_merge_vcf:
#     input:
#         expand(config["dir_data"] + "{cohort}/pindel/chroms/{cohort}.{ref_name}.{tech}.{chrom}.pindel.SV.raw.vcf",
#             contig=[f"chr{i}" for i in range(1,23)] + ["chrX", "chrY", "chrM"])
#     output:
#         vcfgz=config["dir_variants"] + "{prefix}{ref_name}{suffix}.pindel.raw.vcf.gz",
#     run:
#         bcftools = config["software"]["bcftools"]
#         shell("{bcftools} concat -Oz -o {output} {input}")


# with open(f"{output.sample}","w") as file:
#     for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
#         file.write(f"{sample}\n")
# # return samples
# bcftools = config["software"]["bcftools"]
# varscan = config["software"]["varscan"]
# shell("echo '{varscan} mpileup2cns {input.mpileup} {params.extra} --min-coverage 6 --min-reads2 3 --variants 1 "
#       "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {output.sample} | "
#       "{bcftools} view -Oz -o {output.vcf} ' >{log} ")
# shell("{varscan} mpileup2cns {input.mpileup} {params.extra} --min-coverage 6 --min-reads2 3 --variants	1 "
#       "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {output.sample} | "
#       "{bcftools} reheader -f {input.ref_idx} --threads {threads}| {bcftools} view  --threads {threads} -Oz -o {output.vcf} 1>>{log} 2>>{log} ")
