# ======================================================================================================================
# Project: variant_calling_cohort
# Script : pbsv.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: pbsv call
# ======================================================================================================================
rule hificnv:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"]
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        config["dir_data"] + "variants_raw/{cohort}/hificnv/samples/{cohort}.{sample}.{ref_name}.{tech}.hificnv.CNV.ok",
    log:
        config["dir_data"] + "variants_raw/{cohort}/hificnv/logs/{cohort}.{sample}.{ref_name}.{tech}",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/hificnv/logs/{cohort}.{sample}.{ref_name}.{tech}",
    threads: get_run_threads("hificnv")
    run:
        ouput_prefix = f"{output}"[:-3]
        hificnv = config["software"]["hificnv"]
        shell("{hificnv}  --bam {input.bam} --ref {input.ref} --exclude {input.exclude} --threads {threads} --output-prefix {ouput_prefix} 2>>{log} 1>>{log} ")
        shell("touch {output}")

