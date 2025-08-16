# ======================================================================================================================
# Project: variant_calling_cohort
# Script : longshot.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: longshot call
# ======================================================================================================================
rule longshot_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/longshot/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.longshot.SNVIndel.raw.vcf.gz",
        vcf=config["dir_data"] + "variants_raw/{cohort}/longshot/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.longshot.SNVIndel.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/longshot/logs/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.SNVIndel.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/longshot/logs/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.SNVIndel.SV.rtime.tsv",
    threads: get_run_threads("longshot_call")
    run:
        longshot = config["software"]["longshot"]
        bcftools = config["software"]["bcftools"]
        workdir = f"{output.vcf}"[:-4]
        if os.path.exists(f"{workdir}"):
            shell("rm -rf {workdir}")
        shell("mkdir -p {workdir}")
        shell("date > {log}")
        shell("cd {workdir} &&  {longshot} -F -r {wildcards.chrom} --bam {input.bam} --ref {input.ref} --out  {output.vcf} 2>>{log} 1>>{log}")
        shell("{bcftools} sort -Oz -o {output.vcfgz} {output.vcf} 2>>{log} 1>>{log}")
