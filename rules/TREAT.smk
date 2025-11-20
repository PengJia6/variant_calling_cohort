# ======================================================================================================================
# Project: variant_calling_cohort
# Script : trgt.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: trgt call
# ======================================================================================================================
rule TREAT_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
        bed=lambda wildcards: config["refs"][wildcards.ref_name]["tandem_repeat"]["TREAT"]
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        config["dir_data"] + "variants_raw/{cohort}/TREAT/samples/{cohort}.{sample}.{ref_name}.{tech}.TREAT.TR.raw.vcf.gz",
    # vcf=config["dir_data"] + "variants_raw/{cohort}/trgt/samples/{cohort}.{sample}.{ref_name}.{tech}.trgt.TR.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/TREAT/logs/{cohort}.{sample}.{ref_name}.{tech}.TREAT.TR.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/TREAT/logs/{cohort}.{sample}.{ref_name}.{tech}.TREAT.TR.rtime.tsv",
    threads: get_run_threads("TREAT_call")
    run:
        workdir = str(output).rstrip(".vcf.gz") + "_tmp"
        shell("mkdir -p {workdir}")
        treat_path = config["software"]["TREAT"]
        otter_path = config["software"]["otter"]
        treat_env = config["software"]["TREAT_env"]
        treat_dir = "/".join(treat_path.split("/")[:-1])
        otter_dir = "/".join(otter_path.split("/")[:-1])
        # pref = f"{workdir}/" + f'{output}'.split("/")[-1].rstrip(".vcf.gz")
        shell("export PATH={treat_dir}:$PATH && "
              " export PATH={treat_env}:$PATH && "
              " export PATH={otter_dir}:$PATH && "
              " export PATH={otter_dir}:$PATH && "
              "cd {workdir} &&"
              "mkdir -p out && rm -rf pp out &&"
              "{treat_path} assembly -b {input.bed} -i {input.bam} -r {input.ref}  -o out -t {threads} -omc 500 2>{log} 1>{log}"
        )
        shell("ln {workdir}/out/sample.vcf.gz {output}")
#

