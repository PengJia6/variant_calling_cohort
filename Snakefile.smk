# ======================================================================================================================
# Project: alignmentSeq
# Script : Snakefile.smk
# Author : Peng Jia
# Date   :  2024/11/7
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: TODO
# ======================================================================================================================


configfile: "confs/variant_config.yaml"  # config
# configfile: "confs/alignment_reference.yaml"  # config
# configfile: "confs/alignment_samples_align.yaml"  # config
# configfile: "confs/alignment_software.yaml"  # config

include: "rules/common.smk"
# samples_info = config["samples_info"]
wildcard_constraints:
    ref_name=f"{config['ref_id']}",
    # aligner="|".join([aligner for tech, info in config["paras"].items() for aligner in info["aligners"]]),
    tech="|".join(config["name_space"]["tech_names"]),
    cohort="|".join(config["name_space"]["cohort_names"]),
    sample="|".join(config["name_space"]["sample_names"]),
    caller="|".join(config["name_space"]["caller_names"]),
# lib="|".join(config["name_space"]["lib_names"]),
# sample_type="|".join(config["name_space"]["sample_types"]),
# contig="|".join(set([i for i in config["refs"][config['ref_id']]["avaliable"]]))

# print()
# print(config["software"])

targets = []
ref_name = config["ref_id"]
# for aligner in config["aligners"]:
for cohort, cohort_info in config["samples"].items():
    samples_info = {}
    for tech, tech_info in cohort_info["path"].items():
        callers = cohort_info["tech"][tech]["caller"]
        samples = [i for i in tech_info]
        for caller in callers:
            if caller in ["varscan", "bcftools", "GATK_HC", "deepvariant"]:
                targets.append(config["dir_data"] + f"variants_raw/{cohort}/{caller}/{cohort}.{ref_name}.{tech}.{caller}.SNVIndel.raw.vcf.gz.tbi")
            if caller in ["pindel", "manta"]:
                targets.append(config["dir_data"] + f"variants_raw/{cohort}/{caller}/{cohort}.{ref_name}.{tech}.{caller}.SV.raw.vcf.gz.tbi")
            if caller in ["sniffles", "sniffles2", "cutesv", "svision", "svision_pro", "debreak", "pbsv"]:
                for sample in samples:
                    # "samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.vcf.gz","
                    targets.append(config["dir_data"] + f"variants_raw/{cohort}/{caller}/samples/{cohort}.{sample}.{ref_name}.{tech}.{caller}.SV.raw.vcf.gz.tbi")
            if caller in ["hificnv"]:
                for sample in samples:
                    targets.append(config["dir_data"] + f"variants_raw/{cohort}/{caller}/samples/{cohort}.{sample}.{ref_name}.{tech}.{caller}.CNV.ok")
            if caller in ["trgt", ]:
                targets.append(config["dir_data"] + f"variants_raw/{cohort}/trgt/{cohort}.{ref_name}.{tech}.trgt.TR.raw.vcf.gz.tbi",)
            if caller in ["longtr"]:
                for sample in samples:
                    targets.append(config["dir_data"] + f"variants_raw/{cohort}/{caller}/samples/{cohort}.{sample}.{ref_name}.{tech}.{caller}.TR.raw.vcf.gz.tbi")

include: "rules/varscan.smk"
include: "rules/bcftools.smk"
include: "rules/gatk_hc.smk"
include: "rules/deepvariant.smk"
include: "rules/pindel.smk"
include: "rules/sniffles.smk"
include: "rules/cutesv.smk"
include: "rules/svision.smk"
include: "rules/debreak.smk"
include: "rules/pbsv.smk"
include: "rules/longtr.smk"
include: "rules/trgt.smk"
include: "rules/hificnv.smk"

# include: "rules/bam_merge.smk"
# include: "rules/leftalign.smk"
# include: "rules/markdup.smk"
# include: "rules/bwa.smk"
# include: "rules/bowtie2.smk"
# include: "rules/minimap2.smk"
# include: "rules/winnowmap2.smk"
# include: "rules/ngmlr.smk"
# include: "rules/pbmm2.smk"
# include: "rules/lra.smk"

# localrules: all

rule all:
    input:
        targets
