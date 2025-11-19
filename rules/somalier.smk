# ======================================================================================================================
# Project: variant_calling_cohort
# Script : svision.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: svision call
# ======================================================================================================================
rule somalier_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
        sites=lambda wildcards: config["refs"][wildcards.ref_name]["somalier"]
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        config["dir_data"] + "variants_raw/{cohort}/somalier/samples/{cohort}.{sample}.{ref_name}.{tech}.somalier",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/somalier/logs/{cohort}.{sample}.{ref_name}.{tech}.somalier.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/somalier/logs/{cohort}.{sample}.{ref_name}.{tech}.somalier.rtime.tsv",
    threads: get_run_threads("svim_call")
    run:
        somalier = config["software"]["somalier"]
        shell("date > {log}")
        work_dir = "/".join(f"{output}".split("/")[:-1])
        work_dir_tmp = f"{output}_tmp"
        prefix = f"{output}".split("/")[-1][:-9]
        # print(work_dir)
        # print(prefix)
        shell("mkdir -p {work_dir_tmp}")
        shell("cd {work_dir_tmp} && {somalier} extract -d {work_dir_tmp}/ --sites {input.sites} -f  {input.ref} {input.bam}  2>>{log} 1>>{log}")
        shell("ln {work_dir_tmp}/*somalier {output}")
        # shell("rm -rf {work_dir_tmp}")


def get_cohort_samples4somalier(wildcards):
    # for cohort, cohort_info in config["samples"].items():
    #     samples_info = {}
    sample_sites_info = []
    for tech, tech_info in config["samples"][wildcards.cohort]["path"].items():
        samples = [i for i in tech_info]
        for sample in samples:
            sample_sites_info.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/somalier/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{tech}.somalier")

    return sample_sites_info


rule somalier_real:
    input:
        get_cohort_samples4somalier
    output:
        config["dir_data"] + "variants_raw/{cohort}/somalier/{cohort}.{ref_name}.somalier.ok"
    run:
        somalier = config["software"]["somalier"]
        prefix = f"{output}"[:-9]
        shell("{somalier} relate -o {prefix}  {input} ")
        shell("touch {output}")

