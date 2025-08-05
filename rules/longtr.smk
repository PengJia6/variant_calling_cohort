# ======================================================================================================================
# Project: variant_calling_cohort
# Script : longtr.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: longtr call
# ======================================================================================================================
rule longtr_call:
    input:
        unpack(get_samples_bam),
        # bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
        bed=lambda wildcards: config["refs"][wildcards.ref_name]["tandem_repeat"]["LongTR"]
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        # config["dir_data"] + "variants_raw/{cohort}/longtr/samples/chroms/{cohort}.{sample}.{ref_name}.{tech}.longtr.{chrom}.TR.raw.vcf.gz",
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/longtr/chroms/{cohort}.{ref_name}.{tech}.longtr.{chrom}.TR.raw.vcf.gz"
        # vcfgz=config["dir_data"] + "variants_raw/{cohort}/hipstr/chroms/{cohort}.{ref_name}.{tech}.hipstr.{chrom}.TR.raw.vcf.gz",
    log:
        config["dir_data"] + "variants_raw/{cohort}/longtr/logs/{cohort}.{ref_name}.{tech}.{chrom}.longtr.TR.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/longtr/logs/{cohort}.{ref_name}.{tech}.{chrom}.longtr.TR.rtime.tsv",
    threads: get_run_threads("longtr_call")
    run:
        # max_tr_len = config["max_tr_len"]
        bam_str = ",".join([i for i in input.bams])

        longtr = config["software"]["longtr"]
        bcftools = config["software"]["bcftools"]
        vcf_tmp = f"{output.vcfgz}"[:-6] + "tmp.vcf.gz"
        shell("{longtr} --bams {bam_str} --fasta {input.ref} --regions {input.bed} --tr-vcf {vcf_tmp} --chrom {wildcards.chrom}	--output-gls "
              "--max-tr-len 10000 --log {log}")
        shell("{bcftools} sort  -o {output.vcfgz} {vcf_tmp}")


# workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
# bcftools = config["software"]["bcftools"]
# minimap2 = config["software"]["minimap2"]
# wtdbg2 = config["software"]["wtdbg2"]
# prefix = ""
# for i in [longtr, minimap2, wtdbg2]:
#     prefix += "/".join(i.split("/")[:-1])
# if os.path.exists(f"{workdir}"):
#     shell("rm -rf {workdir}")
# shell("mkdir -p {workdir}")
# shell("date > {log}")
# shell("""
#      echo  {longtr} -t {threads} --bam {input.bam} -o {workdir}/ --rescue_large_ins --rescue_dup --poa --ref {input.ref} --min_support 3 2>>{log} 1>>{log}
#         """)
# shell("""
# export PATH={prefix}$PATH &&
# {longtr} -t {threads} --bam {input.bam} -o {workdir}/ --rescue_large_ins --rescue_dup --poa --ref {input.ref} --min_support 3 2>>{log} 1>>{log}
# """)
# shell("date > {log}")
# shell("{bcftools} sort  -o {output.vcf} {workdir}/*vcf")
# shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")
