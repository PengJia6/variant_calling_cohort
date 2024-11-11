# ======================================================================================================================
# Project: alignmentSeq
# Script : Snakefile.smk
# Author : Peng Jia
# Date   :  2024/11/7
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: TODO
# ======================================================================================================================


configfile: "confs/alignment_config.yaml"  # config
# configfile: "confs/alignment_reference.yaml"  # config
# configfile: "confs/alignment_samples_align.yaml"  # config
# configfile: "confs/alignment_software.yaml"  # config

include: "rules/common.smk"
# samples_info = config["samples_info"]
wildcard_constraints:
    ref_name=f"{config['ref_id']}",
    aligner="|".join([aligner for tech, info in config["paras"].items() for aligner in info["aligners"]]),
    tech="|".join(config["name_space"]["tech_names"]),
    cohort="|".join(config["name_space"]["cohort_names"]),
    sample="|".join(config["name_space"]["sample_names"]),
    lib="|".join(config["name_space"]["lib_names"]),
    sample_type="|".join(config["name_space"]["sample_types"]),
    contig="|".join(set([i for i in config["refs"][config['ref_id']]["avaliable"]]))

# print()
print(config["software"])

targets = []

# for aligner in config["aligners"]:
for cohort, cohort_info in config["samples"].items():
    for sample, info in cohort_info.items():
        for sample_type, info_t in info["path"].items():
            for tech in info_t:
                # if tech in config["tech_short_reads"]:
                aligners = config["paras"][tech]["aligners"]
                bam_prefix = "{tech}.{ref_name}.{aligner}"
                if ("mark_duplication" in config["paras"][tech]) and config["paras"][tech]["mark_duplication"]:
                    bam_prefix = bam_prefix + "." + "markdup"
                if ("bam_left_align" in config["paras"][tech] and config["paras"][tech]["bam_left_align"]):
                    bam_prefix = bam_prefix + "." + "leftAlign"
                for aligner in aligners:
                    this_prefix = bam_prefix.format(aligner=aligner,tech=tech,ref_name=config['ref_id'])
                    # config["dir_data"] + "{cohort}/{sample}/temp/{sample}.{sample_type}.{tech}.{ref_name}/{cohort}.{sample}." + \
                    # "{sample_type}_{subsample}.{ref_name}.{tech}.bwa.sorted.addrg.bam",
                    targets.append(config["dir_data"] + f"{cohort}/{sample}/{cohort}.{sample}.{sample_type}.{this_prefix}.bam.bai")

include: "rules/common.smk"
include: "rules/bam_merge.smk"
include: "rules/leftalign.smk"
include: "rules/markdup.smk"
include: "rules/bwa.smk"
include: "rules/bowtie2.smk"
include: "rules/minimap2.smk"
include: "rules/winnowmap2.smk"
include: "rules/ngmlr.smk"
include: "rules/pbmm2.smk"
include: "rules/lra.smk"

# localrules: all

rule all:
    input:
        targets
