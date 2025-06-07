# ======================================================================================================================
# Project: variant_calling_cohort
# Script : trgt.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: trgt call
# ======================================================================================================================
rule trgt_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
        bed=lambda wildcards: config["refs"][wildcards.ref_name]["tandem_repeat"]["TRGT"]
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/trgt/samples/{cohort}.{sample}.{ref_name}.{tech}.trgt.TR.raw.vcf.gz",
    # vcf=config["dir_data"] + "variants_raw/{cohort}/trgt/samples/{cohort}.{sample}.{ref_name}.{tech}.trgt.TR.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/trgt/logs/{cohort}.{sample}.{ref_name}.{tech}.trgt.TR.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/trgt/logs/{cohort}.{sample}.{ref_name}.{tech}.trgt.TR.rtime.tsv",
    threads: get_run_threads("trgt_call")
    run:
        workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
        shell("mkdir -p {workdir}")
        trgt = config["software"]["trgt"]
        bcftools = config["software"]["bcftools"]
        pref = f"{workdir}/" + f'{output}'.split("/")[-1].rstrip(".vcf.gz")
        shell("{trgt} genotype --genome {input.ref} --repeats {input.bed} --reads {input.bam} --output-prefix {pref} --threads {threads} 2>{log} 1>{log}")
        # shell("ln {pref}.vcf.gz {output.vcfgz}")
        shell("{bcftools} sort  -o {output.vcfgz} {pref}.vcf.gz")

        # shell("touch {output}")


def get_trgt_merge_input(wildcards):
    vcfs = []
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        # bams.append(path_bam)
        # bais.append(path_bam + ".bai")
        vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/trgt/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.trgt.TR.raw.vcf.gz")
    return vcfs


def get_trgt_merge_input_tbi(wildcards):
    vcfs = []
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        # bams.append(path_bam)
        # bais.append(path_bam + ".bai")
        vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/trgt/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.trgt.TR.raw.vcf.gz.tbi")
    return vcfs


rule trgt_merge:
    input:
        vcfs=get_trgt_merge_input,
        vcfs_idx=get_trgt_merge_input_tbi,
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        config["dir_data"] + "variants_raw/{cohort}/trgt/{cohort}.{ref_name}.{tech}.trgt.TR.raw.vcf.gz",
    threads: get_run_threads("__default__")
    run:
        trgt = config["software"]["trgt"]
        # pref = f"{workdir}/" + f'{output}'.split("/")[-1].rstrip(".vcf.gz")
        shell("{trgt} merge --genome {input.ref} -Oz -o {output} --vcf {input.vcfs} --force-single ")
# shell("touch {output}")
