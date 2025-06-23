# ======================================================================================================================
# Project: variant_calling_cohort
# Script : varscan.smk TODO check 
# Author : Peng Jia
# Date   :  2024/11/18
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: TODO
# ======================================================================================================================

def get_samples_bam(wildcards):
    bams = []
    bais = []
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        bams.append(path_bam)
        bais.append(path_bam + ".bai")
    return {"bams": bams, "bais": bais,
            "ref": config["refs"][wildcards.ref_name]["fasta"],
            "ref_idx": config["refs"][wildcards.ref_name]["fasta"] + ".fai"
            }


# def get_sample_names(wildcards):


rule samtools_mpileup:
    input:
        unpack(get_samples_bam)
    output:
        config["dir_data"] + "variants_raw/{cohort}/samtools/chroms/{cohort}.{ref_name}.{tech}.{chrom}.samtools.mpileup"
    params:
        extra="",
        dp=5
    log:
        config["dir_data"] + "variants_raw/{cohort}/samtools/logs/{cohort}.{ref_name}.{tech}.{chrom}.samtools.mpileup.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/samtools/logs/{cohort}.{ref_name}.{tech}.{chrom}.samtools.mpileup.rtime.tsv"
    threads: get_run_threads("__default__")
    run:



        samtools = config["software"]["samtools"]
        shell("echo {samtools} mpileup -r {wildcards.chrom} -B -f {input.ref} -o {output} {input.bams} "
              ">{log} ")
        shell("{samtools} mpileup -r {wildcards.chrom} -B -f {input.ref} -o {output} {input.bams} "
              "2>>{log} 1>>{log} ")

#
rule varscan_call_snp_indel:
    input:
        mpileup=config["dir_data"] + "variants_raw/{cohort}/samtools/chroms/{cohort}.{ref_name}.{tech}.{chrom}.samtools.mpileup",
        ref_idx=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai"
    output:
        vcf=config["dir_data"] + "variants_raw/{cohort}/varscan/chroms/{cohort}.{ref_name}.{tech}.varscan.{chrom}.SNVIndel.raw.vcf.gz",
        sample=config["dir_data"] + "variants_raw/{cohort}/varscan/chroms/{cohort}.{ref_name}.{tech}.varscan.{chrom}.SNVIndel.samples.txt"
    params:
        extra="",
    # samples=get_sample_names
    threads: get_run_threads("varscan_call_snp_indel")
    log:
        config["dir_data"] + "variants_raw/{cohort}/varscan/logs/{cohort}.{ref_name}.{tech}.varscan.{chrom}.SNVIndel.varscan.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/varscan/logs/{cohort}.{ref_name}.{tech}.varscan.{chrom}.SNVIndel.varscan.rtime.tsv"
    run:
        with open(f"{output.sample}","w") as file:
            for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
                file.write(f"{sample}\n")
        bcftools = config["software"]["bcftools"]
        varscan = config["software"]["varscan"]
        shell("echo '{varscan} mpileup2cns {input.mpileup} {params.extra} --min-coverage 6 --min-reads2 3 --variants 1 "
              "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {output.sample} | "
              "{bcftools} view -Oz -o {output.vcf} ' >{log} ")
        shell("{varscan} mpileup2cns {input.mpileup} {params.extra} --min-coverage 6 --min-reads2 3 --variants	1 "
              "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {output.sample} | "
              "{bcftools} reheader -f {input.ref_idx} --threads {threads}| {bcftools} view  --threads {threads} -Oz -o {output.vcf} 1>>{log} 2>>{log} ")
