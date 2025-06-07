# ======================================================================================================================
# Project: variant_calling_cohort
# Script : sniffles.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: sniffles and sniffles2 call
# ======================================================================================================================
rule sniffles_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.vcf.gz",
        vcf=config["dir_data"] + "variants_raw/{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.vcf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/sniffles/logs/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/sniffles/logs/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.rtime.tsv",
    threads: get_run_threads("sniffles_call")
    run:
        workdir = str(output.vcf)[:-7] + "_tmp"
        # vcf = f"{output.vcfg}"[:-3]
        # shell("mkdir -p {workdir}")
        vcf_tmp=f"{output.vcf}.tmp.vcf"
        vcf_tmp2=f"{output.vcf}.tmp2.vcf"
        sniffles = config["software"]["sniffles"]
        bcftools = config["software"]["bcftools"]
        shell("date > {log}")
        shell("echo {sniffles} -s 3 -m -i {input.bam} -v {vcf_tmp} -n 3 --cluster -f 0.2 --min_homo_af 0.75 --min_het_af 0.2 -t {threads} 2>>{log} 1>>{log}")
        shell("{sniffles} -s 3  -m {input.bam} -v {vcf_tmp} -n 3 -f 0.2 --min_homo_af 0.75 --min_het_af 0.2 -t {threads} 2>>{log} 1>>{log}")
        # shell("{bcftools} sort -Oz -o {output.vcfgz} {output.vcf}  && rm -rf {workdir}")
        shell("""
        {bcftools} annotate --header-line '##FILTER=<ID=STRANDBIAS,Description="Strand bias detected">' -o {vcf_tmp2} {vcf_tmp}
        """)
        shell("{bcftools} sort  -o {output.vcf} {vcf_tmp2}")
        shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")


rule sniffles2_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/sniffles2/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles2.SV.raw.vcf.gz",
        vcf=config["dir_data"] + "variants_raw/{cohort}/sniffles2/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles2.SV.raw.vcf",
        snf=config["dir_data"] + "variants_raw/{cohort}/sniffles2/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles2.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/sniffles2/logs/{cohort}.{sample}.{ref_name}.{tech}.sniffles2.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/sniffles2/logs/{cohort}.{sample}.{ref_name}.{tech}.sniffles2.SV.rtime.tsv",
    threads: get_run_threads("sniffles2_call")
    run:
        sniffles2 = config["software"]["sniffles2"]
        bcftools = config["software"]["bcftools"]
        vcf_tmp=f"{output.vcf}.tmp.vcf"

        shell("date >{log}")
        shell("echo {sniffles2} --minsupport 3 -i {input.bam} -v {output.vcf} --snf {output.snf} --reference {input.ref} "
              " -t {threads} --minsvlen 10 --mapq 20  --qc-coverage 5   2>>{log} 1>>{log}")
        shell("{sniffles2} --minsupport 3 -i {input.bam} -v {vcf_tmp} --snf {output.snf} --reference {input.ref} "
              " -t {threads} --minsvlen 10 --mapq 20  --qc-coverage 5   2>>{log} 1>>{log}")
        shell("{bcftools} sort  -o {output.vcf} {vcf_tmp}")
        shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")
        # shell("{bcftools} sort -Oz -o {output.vcfgz} {output.vcf}")