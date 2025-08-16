# ======================================================================================================================
# Project: variant_calling_cohort
# Script : clair3.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: clair3 call
# ======================================================================================================================
rule clair3_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/clair3/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.clair3.SNVIndel.raw.vcf.gz",
        # vcf=config["dir_data"] + "variants_raw/{cohort}/clair3/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.clair3.SNVIndel.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/clair3/logs/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.SNVIndel.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/clair3/logs/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.SNVIndel.SV.rtime.tsv",
    threads: get_run_threads("clair3_call")
    run:
        clair3 = config["software"]["clair3"]
        clair3_prefix = "/".join(f"{clair3}".split("/")[:-1])
        bcftools = config["software"]["bcftools"]
        if f"{wildcards.tech}".upper() in ["ILM", "NGS", "BGI", "MGI"]:
            model_path = f"{clair3_prefix}/models/ilmn"
            platform = "ilmn"
        elif f"{wildcards.tech}".upper() in ["ONT", "ONTUL", "ULONT"]:
            model_path = f"{clair3_prefix}/models/ont"
            platform = "ont"
        elif f"{wildcards.tech}".upper() in ["HIFI", "CCS", ]:
            model_path = f"{clair3_prefix}/models/hifi"
            platform = "hifi"
        else:
            print("NO available model provided")
            exit()
        workdir = f"{output.vcfgz}"[:-7]

        if os.path.exists(f"{workdir}"):
            shell("rm -rf {workdir}")
        shell("mkdir -p {workdir}")
        shell("date > {log}")
        shell("cd {workdir} && export PATH={clair3_prefix}:$PATH &&"
              " {clair3} -b {input.bam} -f {input.ref} -m {model_path} -t {threads} -p {platform} -o {workdir} --gvcf --ctg_name={wildcards.chrom} "
              " 2>>{log} 1>>{log}")
        shell("{bcftools} sort -Oz -o {output.vcfgz} {workdir}/merge_output.vcf.gz 2>>{log} 1>>{log}")
