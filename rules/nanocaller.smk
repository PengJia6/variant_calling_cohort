# ======================================================================================================================
# Project: variant_calling_cohort
# Script : nanocaller.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: longshot call
# ======================================================================================================================
rule nanocaller_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/nanocaller/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.nanocaller.SNVIndel.raw.vcf.gz",
        # vcf=config["dir_data"] + "variants_raw/{cohort}/nanocaller/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.nanocaller.SNVIndel.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/nanocaller/logs/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.nanocaller.SNVIndel.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/nanocaller/logs/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.nanocaller.SNVIndel.SV.rtime.tsv",
    threads: get_run_threads("nanocaller_call")
    run:
        # if f"{wildcards.tech}".upper() in ["ILM", "NGS", "BGI", "MGI"]:
        #     # model_path = f"{clair3_prefix}/models/ilmn"
        #     platform = "ilmn"
        if f"{wildcards.tech}".upper() in ["ONT", "ONTUL", "ULONT"]:
            # model_path = f"{clair3_prefix}/models/ont"
            platform = "ont"
            preset = "ont"
        elif f"{wildcards.tech}".upper() in ["HIFI", "CCS", ]:
            # model_path = f"{clair3_prefix}/models/hifi"
            platform = "pacbio"
            preset = "ccs"
        elif f"{wildcards.tech}".upper() in ["CLR", ]:
            # model_path = f"{clair3_prefix}/models/hifi"
            platform = "pacbio"
            preset = "clr"

        else:
            print("NO available model provided")
            exit()
        nanocaller = config["software"]["nanocaller"]
        nanocaller_prefix = "/".join(nanocaller.split("/")[:-1])
        bcftools = config["software"]["bcftools"]
        workdir = f"{output.vcfgz}"[:-7]
        if os.path.exists(f"{workdir}"):
            shell("rm -rf {workdir}")
        shell("mkdir -p {workdir}")
        shell("date > {log}")
        shell("cd {workdir} && "
              "export PATH={nanocaller_prefix}:$PATH && "
              "{nanocaller} --regions {wildcards.chrom} --cpu {threads}  --sequencing {platform} --mode all "
              "--preset {preset} --bam {input.bam} --ref {input.ref} --sample {wildcards.cohort}_{wildcards.sample}   "
              "2>>{log} 1>>{log}")
        # shell ("ln -s {workdir}/variant_calls.vcf.gz {output.vcf}")
        shell("{bcftools} sort -Oz -o {output.vcfgz} {workdir}/variant_calls.vcf.gz 2>>{log} 1>>{log}")
