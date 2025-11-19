# ======================================================================================================================
# Project: variant_calling_cohort
# Script : trgt.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: trgt call
# ======================================================================================================================
rule EH_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
        bed=lambda wildcards: config["refs"][wildcards.ref_name]["tandem_repeat"]["ExpansionHunter"]
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/ExpansionHunter/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.ExpansionHunter.TR.raw.vcf.gz",
        # vcfgz=config["dir_data"] + "variants_raw/{cohort}/gangstr/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.gangstr.TR.raw.vcf.gz",

    # vcf=config["dir_data"] + "variants_raw/{cohort}/trgt/samples/{cohort}.{sample}.{ref_name}.{tech}.trgt.TR.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/ExpansionHunter/logs/{cohort}.{sample}.{ref_name}.{tech}.{chrom}/ExpansionHunter.TR.log",
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/ExpansionHunter/logs/{cohort}.{sample}.{ref_name}.{tech}.{chrom}/ExpansionHunter.TR.rtime.tsv",
    threads: get_run_threads("ExpansionHunter_call")
    run:
        workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
        shell("mkdir -p {workdir}")
        this_bed=f"{input.bed}"[:-3]+f"/{wildcards.chrom}.EH.json"
        ExpansionHunter = config["software"]["ExpansionHunter"]
        bcftools = config["software"]["bcftools"]
        pref = f"{workdir}/" + f'{output}'.split("/")[-1].rstrip("vcf.gz")
        shell("cd {workdir} && {ExpansionHunter} --reads {input.bam} --reference {input.ref} "
              "--variant-catalog {this_bed} "
              "--output-prefix {pref} -n {threads} 2>{log} 1>{log}")
        # shell("ln {pref}.vcf.gz {output.vcfgz}")
        shell("{bcftools} view -Oz  -o {output.vcfgz} {pref}.vcf")


def get_eh_merge_input(wildcards):
    vcfs = []
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/ExpansionHunter/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.ExpansionHunter.TR.raw.vcf.gz")
    return vcfs


def get_eh_merge_input_tbi(wildcards):
    vcfs = []
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        # bams.append(path_bam)
        # bais.append(path_bam + ".bai")
        vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/ExpansionHunter/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.ExpansionHunter.TR.raw.vcf.gz.tbi")
    return vcfs


rule EH_merge:
    input:
        vcfs=get_eh_merge_input,
        vcfs_idx=get_eh_merge_input_tbi,
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        config["dir_data"] + "variants_raw/{cohort}/ExpansionHunter/{cohort}.{ref_name}.{tech}.ExpansionHunter.TR.raw.vcf.gz",
    threads: get_run_threads("__default__")
    run:
        trgt = config["software"]["trgt"]
        # pref = f"{workdir}/" + f'{output}'.split("/")[-1].rstrip(".vcf.gz")
        shell("{trgt} merge --genome {input.ref} -Oz -o {output} --vcf {input.vcfs} --force-single ")
# shell("touch {output}")
