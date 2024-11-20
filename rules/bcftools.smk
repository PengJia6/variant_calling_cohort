# ======================================================================================================================
# Project: variant_calling_cohort
# Script : varscan.smk TODO check 
# Author : Peng Jia
# Date   :  2024/11/18
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: TODO
# ======================================================================================================================



rule bcftools_call:
    input:
        unpack(get_samples_bam)
    output:
        vcf = config["dir_data"] + "{cohort}/bcftools/chroms/{cohort}.{ref_name}.{tech}.bcftools.{chrom}.SNVIndel.raw.vcf.gz",
    params:
        extra="",
        dp=5
    log:
        config["dir_data"] + "{cohort}/bcftools/logs/{cohort}.{ref_name}.{tech}.{chrom}.bcftools.mpileup.log"
    benchmark:
        config["dir_data"] + "{cohort}/bcftools/logs/{cohort}.{ref_name}.{tech}.{chrom}.bcftools.mpileup.rtime.tsv"
    threads: get_run_threads("bcftools_call")
    run:
        bcftools = config["software"]["bcftools"]
        # shell("echo {bcftools} mpileup -r {wildcards.chrom}  --threads {threads} -B -f {input.ref} -o {output} {input.bams} "
        #       ">{log} ")
        shell("{bcftools} mpileup -r {wildcards.chrom} --threads {threads} -f {input.ref} {input.bams} |"
              "{bcftools} call --threads {threads} -mv -Oz -o {output} "
              "2>>{log} 1>>{log} ")


#
# rule bcftools_call:
#     input:
#         config["dir_data"] + "{cohort}/samtools/chroms/{cohort}.{ref_name}.{tech}.{chrom}.samtools.mpileup"
#     output:
#         vcf=config["dir_data"] + "{cohort}/bcftools/chroms/{cohort}.{ref_name}.{tech}.bcftools.{chrom}.SNVIndel.raw.vcf.gz",
#         sample=config["dir_data"] + "{cohort}/bcftools/chroms/{cohort}.{ref_name}.{tech}.bcftools.{chrom}.SNVIndel.samples.txt"
#     params:
#         extra="",
#     # samples=get_sample_names
#     threads: get_run_threads("bcftools_call")
#     log:
#         config["dir_data"] + "{cohort}/bcftools/logs/{cohort}.{ref_name}.{tech}.bcftools.{chrom}.SNVIndel.bcftools.log"
#     benchmark:
#         config["dir_data"] + "{cohort}/bcftools/logs/{cohort}.{ref_name}.{tech}.bcftools.{chrom}.SNVIndel.bcftools.rtime.tsv"
#     run:
#         with open(f"{output.sample}","w") as file:
#             for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
#                 file.write(f"{sample}\n")
#         # return samples
#         bcftools = config["software"]["bcftools"]
#         # shell("echo '{varscan} mpileup2cns {input} {params.extra} --min-coverage 6 --min-reads2 3 "
#         #       "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {output.sample} | "
#         #       "{bcftools} view -Oz -o {output.vcf} ' >{log} ")
#         shell("{bcftools} call -mv -Oz -o {output.vcf} --threads {threads} {input}"
#               "{varscan} mpileup2cns {input} {params.extra} --min-coverage 6 --min-reads2 3 "
#               "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {output.sample} 2 >>{log}| "
#               "{bcftools} view -Oz -o {output.vcf} 1>>{log} 2>>{log} ")
#
#
