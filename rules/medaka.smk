# ======================================================================================================================
# Project: variant_calling_cohort
# Script : medaka_var.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: svision call
# ======================================================================================================================
rule medaka_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/medaka/samples/{cohort}.{sample}.{ref_name}.{tech}.medaka.SNVIndel.raw.vcf.gz",
        vcf=config["dir_data"] + "variants_raw/{cohort}/medaka/samples/{cohort}.{sample}.{ref_name}.{tech}.medaka.SNVIndel.raw.vcf",
    # vcf=config["dir_data"] + "variants_raw/{cohort}/nanovar/samples/{cohort}.{sample}.{ref_name}.{tech}.nanovar.SV.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/medaka/logs/{cohort}.{sample}.{ref_name}.{tech}.medaka.SNVIndel.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/medaka/logs/{cohort}.{sample}.{ref_name}.{tech}.medaka.SNVIndel.rtime.tsv",
    threads: get_run_threads("nanovar_call")
    run:
        if f"{wildcards.tech}".upper() in ["ONT", "ONTUL", "ULONT"]:
            # model_path = f"{clair3_prefix}/models/ont"
            preset = "ont"
            workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
            if os.path.exists(f"{workdir}"):
                shell("rm -rf {workdir}")
            shell("mkdir -p {workdir}")

            medaka = config["software"]["medaka"]
            bcftools = config["software"]["bcftools"]
            medaka_prefix = "/".join(medaka.split("/")[:-1])
            # ref_prefix = "/".join(f"{input.ref}".split("/")[:-1])
            # ref_file = f"{input.ref}".split("/")[-1]
            # bgzip = config["software"]["bgzip"]
            shell("date > {log}")
            shell("cd {workdir} && export PATH={medaka_prefix}:$PATH && "
                  "{medaka}_variant -i {input.bam} "
                  "2>>{log} 1>>{log}")
            shell("{bcftools} sort  -o {output.vcf} {workdir}/*.nanovar.pass.vcf")
            shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")
        elif f"{wildcards.tech}".upper() in ["HIFI", "CCS", ]:
            # model_path = f"{clair3_prefix}/models/hifi"
            preset = "ccs"
        elif f"{wildcards.tech}".upper() in ["CLR", ]:
            # model_path = f"{clair3_prefix}/models/hifi"
            preset = "clr"

        else:
            print("NO available model provided")
            exit()

