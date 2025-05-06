# ======================================================================================================================
# Project: variant_calling_cohort
# Script : pbsv.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: pbsv call
# ======================================================================================================================
rule pbSV_discovery:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        sig=config["dir_data"] + "variants_raw/{cohort}/pbsv/samples/{cohort}.{sample}.{ref_name}.{tech}/chroms/{contig}.pbsv.SV.raw.svsig.gz",
    log:
        config["dir_data"] + "variants_raw/{cohort}/pbsv/logs/{cohort}.{sample}.{ref_name}.{tech}/chroms/{contig}.pbsv_dis.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/pbsv/logs/{cohort}.{sample}.{ref_name}.{tech}/chroms/{contig}.pbsv_dis.SV.rtime.tsv",
    threads: get_run_threads("pbsv_call")
    run:
        output_pre = "/".join(str(output).split("/")[:-1])
        shell("mkdir -p {output_pre}")
        pbsv = config["software"]["pbsv"]
        shell("echo {pbsv} discover -s {wildcards.sample} -b {input.ref} --region {wildcards.contig}  {input.bam} {output.sig} 2>>{log} 1>>{log} ")
        shell("{pbsv} discover -s {wildcards.sample} -b {input.ref}  {input.bam} {output.sig} 2>>{log} 1>>{log} ")


def get_pbsv_contig_sig(wildcards):
    sigs = expand(config["dir_data"] + "variants_raw/{cohort}/pbsv/samples/{cohort}.{sample}.{ref_name}.{tech}/chroms/{contig}.pbsv.SV.raw.svsig.gz",
        cohort=wildcards.cohort,ref_name=wildcards.ref_name,sample=wildcards.sample,tech=wildcards.tech,
        contig=config["refs"][wildcards.ref_name]["available_chrom"]
                  )
    return sigs
rule pbsv_call:
    input:
        sigs=get_pbsv_contig_sig,
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/pbsv/samples/{cohort}.{sample}.{ref_name}.{tech}.pbsv.SV.raw.vcf.gz",
        vcf=config["dir_data"] + "variants_raw/{cohort}/pbsv/samples/{cohort}.{sample}.{ref_name}.{tech}.pbsv.SV.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/pbsv/logs/{cohort}.{sample}.{ref_name}.{tech}.pbsv.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/pbsv/logs/{cohort}.{sample}.{ref_name}.{tech}.pbsv.SV.rtime.tsv",
    threads: get_run_threads("__default__")
    run:
        workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
        pbsv = config["software"]["pbsv"]
        bcftools = config["software"]["bcftools"]
        vcf_tmp = f"{output.vcf}.tmp.vcf"
        if os.path.exists(f"{workdir}"):
            shell("rm -rf {workdir}")
        shell("mkdir -p {workdir}")
        shell("date > {log}")
        shell("echo {pbsv} call -t DEL,INS,INV,DUP,BND,CNV --ccs -j {threads} {input.ref} {input.sigs} {vcf_tmp} 2>>{log} 1>>{log}")
        shell("{pbsv} call -t DEL,INS,INV,DUP,BND,CNV --ccs -j {threads} {input.ref} {input.sigs} {vcf_tmp} 2>>{log} 1>>{log}")
        shell("date > {log}")
        shell("{bcftools} sort  -o {output.vcf} {vcf_tmp}")
        shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")
