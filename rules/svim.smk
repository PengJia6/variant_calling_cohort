# ======================================================================================================================
# Project: variant_calling_cohort
# Script : svision.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: svision call
# ======================================================================================================================
rule svim_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/svim/samples/{cohort}.{sample}.{ref_name}.{tech}.svim.SV.raw.vcf.gz",
        vcf=config["dir_data"] + "variants_raw/{cohort}/svim/samples/{cohort}.{sample}.{ref_name}.{tech}.svim.SV.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/svim/logs/{cohort}.{sample}.{ref_name}.{tech}.svim.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/svim/logs/{cohort}.{sample}.{ref_name}.{tech}.svim.SV.rtime.tsv",
    threads: get_run_threads("svim_call")
    run:
        workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
        if os.path.exists(f"{workdir}"):
            shell("rm -rf {workdir}")
        shell("mkdir -p {workdir}")

        svim = config["software"]["svim"]
        bcftools = config["software"]["bcftools"]

        bgzip = config["software"]["bgzip"]
        shell("date > {log}")
        shell("{svim} alignment --min_sv_size 30 --sample {wildcards.sample} {workdir} {input.bam} {input.ref}   2>>{log} 1>>{log}")
        shell("{bcftools} sort  -o {output.vcf} {workdir}/variants.vcf")
        shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")

