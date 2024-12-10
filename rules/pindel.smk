# ======================================================================================================================
# Project: variant_calling_cohort
# Script : pindel.smk
# Author : Peng Jia
# Date   :  2024/11/18
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: TODO
# ======================================================================================================================
import random

import numpy as np
import pysam


def extrac_insert(bam_path):
    sam_file = pysam.Samfile(bam_path)
    # ead_groups = bam_file.header['RG']

    # 遍历每个read group，获取sample信息
    samples = []
    for rg in sam_file.header["RG"]:
        sample = rg.get('SM',None)  # 'SM'是read group中的sample字段
        if sample and sample not in samples:
            samples.append(sample)
    if len(samples) != 1:
        return -1
    insert_size = []
    for chrom, chrom_len in zip(sam_file.header.references,sam_file.header.lengths):
        if chrom_len < 1000000: continue
        print(chrom,chrom_len)
        pos_list = random.sample(range(chrom_len),10)
        for pos in pos_list:
            read_num = 0
            for read in sam_file.fetch(chrom,pos):
                # print(read.mapping_quality)
                if read.is_unmapped or read.is_secondary or read.mapping_quality < 40 or \
                        read.mate_is_unmapped or read.next_reference_name != chrom or \
                        read.is_duplicate:
                    continue
                if read.is_reverse: continue
                if read.mate_is_forward: continue
                insert_size.append(abs(read.reference_start - read.next_reference_start) + read.query_length)
                read_num += 1
                # print(read_num)
                if read_num > 100: break
    sam_file.close()
    return [bam_path, f"{int(np.median(insert_size))}", samples[0]]


rule pindel_call:
    input:
        unpack(get_samples_bam),
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
    output:
        protected(config["dir_data"] + "{cohort}/pindel/chroms/{cohort}.{ref_name}.{tech}.{chrom}.pindel")
    params:
        extra="",
        dp=5
    log:
        config["dir_data"] + "{cohort}/pindel/logs/{cohort}.{ref_name}.{tech}.{chrom}.pindel.mpileup.log"
    benchmark:
        config["dir_data"] + "{cohort}/pindel/logs/{cohort}.{ref_name}.{tech}.{chrom}.pindel.mpileup.rtime.tsv"
    threads: get_run_threads("pindel_call")
    run:
        pindel = config["software"]["pindel"]
        pindel_conf = str(output) + ".conf"
        # sample_name = f"{input.bam}".split("/")[-1][:-4]
        config_file = open(str(output) + ".conf","w")
        for bam_path in input.bams:
            conf_info = extrac_insert(bam_path)
            print(conf_info)
            config_file.write(" ".join(conf_info) + "\n")
        config_file.close()
        pindel_output = str(output).split("/")[-1]
        pindel_prefix = "/".join(str(output).split("/")[:-1])
        pindel = config["software"]["pindel"]
        shell("{pindel} -i {pindel_conf} -f {input.ref} -c {wildcards.chrom} -o {output} -T {threads} -x2 -l -J {input.exclude} 2>{log} 1>{log} ")
        shell("touch {output}")
        #

        # shell("echo {samtools} mpileup -r {wildcards.chrom} -B -f {input.ref} -o {output} {input.bams} "
        #       ">{log} ")
        # shell("{samtools} mpileup -r {wildcards.chrom} -B -f {input.ref} -o {output} {input.bams} "
        #       "2>>{log} 1>>{log} ")

#
rule pindel2vcf:
    input:
        config["dir_data"] + "{cohort}/pindel/chroms/{cohort}.{ref_name}.{tech}.{chrom}.pindel",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        vcf=config["dir_data"] + "{cohort}/pindel/chroms/{cohort}.{ref_name}.{tech}.pindel.{chrom}.SV.raw.vcf"
    # vcf=config["dir_data"] + "{cohort}/varscan/chroms/{cohort}.{ref_name}.{tech}.varscan.{chrom}.SNVIndel.raw.vcf.gz",
    # sample=config["dir_data"] + "{cohort}/varscan/chroms/{cohort}.{ref_name}.{tech}.varscan.{chrom}.SNVIndel.samples.txt"
    params:
        extra="",
    # samples=get_sample_names
    threads: get_run_threads("pindel2vcf")
    log:
        config["dir_data"] + "{cohort}/pindel/logs/{cohort}.{ref_name}.{tech}.pindel.{chrom}.pindel2vcf.log"
    benchmark:
        config["dir_data"] + "{cohort}/pindel/logs/{cohort}.{ref_name}.{tech}.pindel.{chrom}.pidnel2vcf.rtime.tsv"
    run:
        pindel2vcf = config["software"]["pindel2vcf"]
        out_prefix = str(output.vcf)[:-4]
        shell("{pindel2vcf} -P {input} -r {input} -R {wildcards.ref_name} -d 2024 -v {output.vcf} -is 30 -as 100000000 -b -e 10 -ss 5    ")

rule pidnel_reheadervcf:
    input:
        vcf=config["dir_data"] + "{cohort}/pindel/chroms/{cohort}.{ref_name}.{tech}.pindel.{chrom}.SV.raw.vcf",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],

    output:
        vcf=config["dir_data"] + "{cohort}/pindel/chroms/{cohort}.{ref_name}.{tech}.pindel.{chrom}.SV.raw.vcf.gz"
    wildcard_constraints:
        contig="|".join([f"chr{i}" for i in range(1,23)] + ["chrX", "chrY", "chrM"])
    run:
        bcftools = config["software"]["bcftools"]
        shell("{bcftools} reheader --fai {input.ref} -Oz {input.vcf} > {output}")


# def get_vcfs_for_pindel_merge_vcf(wildcards):
#     gvcfs = []
#     vcfs = []
#     # for cohort, cohort_info in.items():
#     for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
#         # for sample_name, info_t in config["samples_info"][f"{wildcards.cohort}"].items():
#         #     for sample_type, info_tt in info_t["path"].items():
#         # print(wildcards.cohort,sample_type,sample_name,wildcards.this_prefix)
#         # for tech, info in info_tt.items():
#         vcfs.append(config["dir_data"] + f"{wildcards.cohort}/deepvariant/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.deepvariant.SNVIndel.vcf.gz")
#         gvcfs.append(config["dir_data"] + f"{wildcards.cohort}/deepvariant/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.deepvariant.SNVIndel.g.vcf.gz")
#     return {"gvcf_gz": gvcfs, "vcf_gz": vcfs}
# rule pindel_merge_vcf:
#     input:
#         expand(config["dir_data"] + "{cohort}/pindel/chroms/{cohort}.{ref_name}.{tech}.{chrom}.pindel.SV.raw.vcf",
#             contig=[f"chr{i}" for i in range(1,23)] + ["chrX", "chrY", "chrM"])
#     output:
#         vcfgz=config["dir_variants"] + "{prefix}{ref_name}{suffix}.pindel.raw.vcf.gz",
#     run:
#         bcftools = config["software"]["bcftools"]
#         shell("{bcftools} concat -Oz -o {output} {input}")


# with open(f"{output.sample}","w") as file:
#     for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
#         file.write(f"{sample}\n")
# # return samples
# bcftools = config["software"]["bcftools"]
# varscan = config["software"]["varscan"]
# shell("echo '{varscan} mpileup2cns {input.mpileup} {params.extra} --min-coverage 6 --min-reads2 3 --variants 1 "
#       "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {output.sample} | "
#       "{bcftools} view -Oz -o {output.vcf} ' >{log} ")
# shell("{varscan} mpileup2cns {input.mpileup} {params.extra} --min-coverage 6 --min-reads2 3 --variants	1 "
#       "--min-avg-qual 8  --min-var-freq 0.1 --output-vcf 1 --vcf-sample-list {output.sample} | "
#       "{bcftools} reheader -f {input.ref_idx} --threads {threads}| {bcftools} view  --threads {threads} -Oz -o {output.vcf} 1>>{log} 2>>{log} ")
