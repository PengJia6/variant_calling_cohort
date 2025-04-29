# ======================================================================================================================
# Project: variant_calling_cohort
# Script : svision.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: svision call
# ======================================================================================================================
rule svision_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/svision/samples/{cohort}.{sample}.{ref_name}.{tech}.svision.SV.raw.vcf.gz",
        vcf=config["dir_data"] + "variants_raw/{cohort}/svision/samples/{cohort}.{sample}.{ref_name}.{tech}.svision.SV.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/svision/logs/{cohort}.{sample}.{ref_name}.{tech}.svision.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/svision/logs/{cohort}.{sample}.{ref_name}.{tech}.svision.SV.rtime.tsv",
    threads: get_run_threads("svision_call")
    run:
        workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
        if os.path.exists(f"{workdir}"):
            shell("rm -rf {workdir}")
        my_model = config["software"]["SVision_model"]
        SVision = config["software"]["SVision"]
        bcftools = config["software"]["bcftools"]

        bgzip = config["software"]["bgzip"]
        shell("date > {log}")
        shell("echo {SVision} -t {threads} -o {workdir} -n {wildcards.sample} -b {input.bam} "
              "-g {input.ref} -m {my_model} --min_sv_size 30 -s10  2>>{log} 1>>{log}")
        shell("{SVision} -t {threads} -o {workdir} -n {wildcards.sample} -b {input.bam} "
              "-g {input.ref} -m {my_model} --min_sv_size 30 -s10  2>>{log} 1>>{log}")
        # shell("{bgzip} -c {workdir}/{wildcards.sample}.svision.s10.vcf > {output.vcfgz}")
        # shell("mv {workdir}/{wildcards.sample}.svision.s10.vcf {output.vcf}")
        shell("{bcftools} sort  -o {output.vcf} {workdir}/*.vcf")
        shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")


rule svision_pro_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/svision_pro/samples/{cohort}.{sample}.{ref_name}.{tech}.svision_pro.SV.raw.vcf.gz",
        vcf=config["dir_data"] + "variants_raw/{cohort}/svision_pro/samples/{cohort}.{sample}.{ref_name}.{tech}.svision_pro.SV.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/svision_pro/logs/{cohort}.{sample}.{ref_name}.{tech}.svision_pro.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/svision_pro/logs/{cohort}.{sample}.{ref_name}.{tech}.svision_pro.SV.rtime.tsv",
    threads: get_run_threads("svision_pro_call")
    run:
        workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
        if os.path.exists(f"{workdir}"):
            shell("rm -rf {workdir}")
        my_model = config["software"]["SVision_pro_model"]
        SVision_pro = config["software"]["SVision_pro"]
        bcftools = config["software"]["bcftools"]

        bgzip = config["software"]["bgzip"]
        if f"{wildcards.tech}".upper() in ["HIFI","CCS"]:
            preset="hifi"
        elif f"{wildcards.tech}".upper() in ["ONT","ONTUL","ULONT"]:
            preset= "error-prone"
        elif f"{wildcards.tech}".upper() in ["ASSM", "ASM"]:
            preset = "asm"

        shell("date > {log}")
        shell("echo {SVision_pro} --preset {preset} --process_num {threads} --target_path {input.bam} --genome_path {input.ref} --model_path {my_model} "
              "--out_path {workdir} --sample_name {wildcards.sample} --detect_mode germline "
              "--min_supp 2 --min_sv_size 30 --max_sv_size 100000 2>> {log} 1>>{log}"
              )
        shell("{SVision_pro} --preset {preset} --process_num {threads} --target_path {input.bam} --genome_path {input.ref} --model_path {my_model} "
              "--out_path {workdir} --sample_name {wildcards.sample} --detect_mode germline "
              "--min_supp 2 --min_sv_size 30 --max_sv_size 100000 2>> {log} 1>>{log}"
        )
        shell("{bcftools} sort  -o {output.vcf} {workdir}/*vcf")
        shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")
        # shell("{SVision} -t {threads} -o {workdir} -n {wildcards.sample} -b {input.bam} "
        #       "-g {input.ref} -m {my_model} --min_sv_size 30 -s10  2>>{log} 1>>{log}")
        # shell("{bgzip} -c {workdir}/{wildcards.sample}.svision.s10.vcf > {output.vcfgz}")
        # shell("mv {workdir}/{wildcards.sample}.svision.s10.vcf {output.vcf}")
        # shell("touch {output.vcf}")

