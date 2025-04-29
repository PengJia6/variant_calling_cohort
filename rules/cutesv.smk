# ======================================================================================================================
# Project: variant_calling_cohort
# Script : cutesv.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: cutesv call
# ======================================================================================================================
rule cutesv_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "{cohort}/cutesv/samples/{cohort}.{sample}.{ref_name}.{tech}.cutesv.SV.raw.vcf.gz",
        vcf=config["dir_data"] + "{cohort}/cutesv/samples/{cohort}.{sample}.{ref_name}.{tech}.cutesv.SV.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "{cohort}/cutesv/logs/{cohort}.{sample}.{ref_name}.{tech}.cutesv.SV.log",

    benchmark:
        config["dir_data"] + "{cohort}/cutesv/logs/{cohort}.{sample}.{ref_name}.{tech}.cutesv.SV.rtime.tsv",
    threads: get_run_threads("cutesv_call")
    run:
        workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
        cuteSV = config["software"]["cuteSV"]
        bcftools = config["software"]["bcftools"]
        if os.path.exists(f"{workdir}"):
            shell("rm -rf {workdir}")
        shell("mkdir -p {workdir}")
        shell("date > {log}")
        shell("echo {cuteSV} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 "
              "-s 2 --genotype "
              "--diff_ratio_merging_DEL	0.5 --max_cluster_bias_DEL	1000 "
              "--threads {threads} --sample	{wildcards.sample} --min_size 20 "
              "{input.bam} {input.ref} {output.vcf} {workdir} 2>>{log} 1>>{log} ")
        shell("{cuteSV} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 "
              "-s 2 --genotype "
              "--diff_ratio_merging_DEL	0.5 --max_cluster_bias_DEL	1000 "
              "--threads {threads} --sample	{wildcards.sample} --min_size 20 "
              "{input.bam} {input.ref} {output.vcf} {workdir} 2>>{log} 1>>{log} ")
        shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf} 2>>{log} 1>>{log}")
