# ======================================================================================================================
# Project: HiFi_analysis_pipeline
# Script : common.py.smk
# Author : Peng Jia
# Date   :  2022/9/29
# Email  : pengjia@stu.xjtu.edu.cn
# Description: common functions for this pipeline
# ======================================================================================================================

import os.path
from shutil import unpack_archive

import yaml


# used by varscan, pindel
def get_samples_bam(wildcards):
    bams = []
    bais = []
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        bams.append(path_bam)
        bais.append(path_bam + ".bai")
    return {"bams": bams, "bais": bais,
            "ref": config["refs"][wildcards.ref_name]["fasta"],
            "ref_idx": config["refs"][wildcards.ref_name]["fasta"] + ".fai"
            }


def get_software(name):
    if "software" not in config:
        print("[Error] no software set in config file.")
    if name in config["software"]:
        return config["software"][name]
    else:
        print(f"[Error] no {name} set in config of software.")


def get_run_threads(rule_name):
    if "threads" not in config:
        print("[Error] no threads set in config file.")
    if rule_name in config["threads"]:
        return config["threads"][rule_name]
    else:
        return config["threads"]["__default__"]


def check_samples(conf_samples):
    error = False
    file_samples = open(conf_samples,"r")
    samples_info = yaml.safe_load(file_samples)
    if "samples" not in config:
        config["samples"] = {}
    for cohort in samples_info:
        config["samples"][cohort] = samples_info[cohort]
    file_samples.close()
    # if "samples" not in config:
    # samples_info = config["samples"]
    # print(samples_info)
    cohort_names = []
    sample_names = []
    # sample_types = []
    tech_names = []
    caller_names = []

    for cohort, cohort_info in samples_info.items():
        cohort_names.append(cohort)
        this_techs = []
        if "path" not in cohort_info or "tech" not in cohort_info:  # "info" not in cohort_info or
            print(f"[Error] The sample config setting of {cohort} is not correct!\n"
                  f"e.g.\n"
                  f"-------------------------\n"
                  f"{cohort}:\n"
                  f"    path:\n"
                  f"        tech1:\n"
                  f"            sample1: path/to/sample1.bam\n"
                  f"    info:\n"
                  f"        sample1:\n"
                  f"            sex:\n"
                  f"            ...\n"
                  f"    tech:\n"
                  f"        tech1:\n"
                  f"            caller:\n"
                  f"               - varscan\n"
                  f"-------------------------\n"

            )
            error = True
        else:
            for tech, info_t in cohort_info["path"].items():
                tech_names.append(tech)
                this_techs.append(tech)
                if tech not in cohort_info["tech"]:
                    print(f"[Error] There is not {tech} info in sample configure file!")
                    error = True
                else:

                    for sample_name, bam_file in info_t.items():
                        sample_names.append(sample_name)
                        if len(bam_file) < 1 or not (isinstance(bam_file,str)) or (not bam_file.endswith(".bam")):
                            print(f"[Error] The bam file for {sample_name} of {tech} is not available!")
                            error = True
                for caller in cohort_info["tech"][tech]["caller"]:
                    if caller not in caller_names:
                        caller_names.append(caller)

    config["name_space"] = {"cohort_names": list(set(cohort_names)),
                            "sample_names": list(set(sample_names)),
                            "caller_names": list(set(caller_names)),
                            "tech_names": list(set(tech_names))}

    return True if not error else False


def check_config_file(config):
    error = False
    if "path_ref_config" not in config:
        config["path_ref_config"] = "confs/alignment_reference.yaml"
    else:
        if not os.path.exists(config["path_ref_config"]):
            print(f"[Error] The config file {config['path_ref_config']} for reference genome "
                  f"('path_ref_config' in config.yaml) is not available! "
                  "Please provide a config file (yaml format) for the reads of the samples in this project! "
                  "You can setting the config file according the template xxxx")  # TODO add the link
            error = True
        else:
            config["refs"] = yaml.safe_load(open(config["path_ref_config"]))["reference"]
    if "ref_id" not in config:
        print(f"[Error] No reference id (ref_id) setting in config file")
        error = True
    else:
        if (config["ref_id"] not in config["refs"]) or ("fasta" not in config["refs"][config["ref_id"]]):
            print(f"[Error] no setting of {config['ref_id']} in {config['path_ref_config']} e.g.\n"
                  f"-----------------------\n"
                  f"{config['ref_id']}:\n\tfasta: /path/to/ref.fasta\n\t...\n"
                  f"-----------------------\n")
            error = True
    if "path_software_config" not in config:
        config["path_software_config"] = "confs/alignment_software.yaml"
    if not os.path.exists(config["path_software_config"]):
        print(f"[Error] The config file {config['path_software_config']} for samples "
              f"('path_software_config' in config.yaml) is not available! "
              "Please provide a config file (yaml format) for the reads of the samples in this project! "
              "You can setting the config file according the template xxxx")  # TODO add the link
        error = True
    else:
        config["software"] = yaml.safe_load(open(config["path_software_config"]))["software"]

    # if "path_cluster_config" not in config:
    #     config["path_cluster_config"] = "confs/alignment_cluster.yaml"
    # if not os.path.exists(config["path_cluster_config"]):
    #     print(f"[Error] The config file {config['path_cluster_config']} for resources"
    #           f" ('path_cluster_config' in config.yaml) is not available! "
    #           "Please provide a config file (yaml format) to ensure that each rule utilizes a sensible threads! "
    #           "You can setting the config file according the template xxxx")  # TODO add the link
    #     error = True
    # else:
    #     config["threads"] = yaml.safe_load(open(config["path_cluster_config"]))

    if "path_samples_config" not in config:
        config["path_samples_config"] = ["confs/alignment_samples.yaml"]
    for sample_config in config["path_samples_config"]:
        if not os.path.exists(sample_config):
            print(
                f"[Error] The config file {sample_config} for samples ('path_samples_config' in config.yaml) is not available! "
                "Please provide a config file (yaml format) for the reads of the samples in this project! "
                "You can setting the config file according the template")
            error = True
        else:
            # check_samples(sample_config)
            # print(samples_info)
            if not check_samples(sample_config):
                error = True
    for soft in config["required_software"]:
        if soft not in config["software"]:
            print(f"[Error] No {soft} setting in {config['path_software_config']}.")
            error = True
    # else:
    # config["samples_info"] = samples_info
    # print()
    # config["samples"] = list([s for c, info in samples_info.items() for s in info])
    # config["sample_types"] = list(
    #     set([i for c, info in samples_info.items() for s, info2 in info.items() for i in info2["path"]]))
    # config["techs"] = list(
    #     set([k for c, info in samples_info.items() for s, info2 in info.items() for i, info3 in
    #          info2["path"].items() for k in info3]))
    # config["cohorts"] = list([c for c in samples_info])

    config["dir_data"] = os.path.abspath("variants_raw") if "dir_data" not in config else config["dir_data"]
    config["dir_data"] = config["dir_data"] if config["dir_data"].endswith("/") else config["dir_data"] + "/"
    # config["dir_ref"] = config["dir_data"] + "ref/"
    # config["dir_data"] = config["dir_data"] + "lr_genome/"
    # config["dir_raw_data"] = config["dir_data"] + "rawdata/"
    # config["dir_aligned_reads"] = config["dir_data"] + "aligned_reads/"
    # config["dir_reports"] = config["dir_data"] + "reports/"
    # config["dir_tmp"] = config["dir_data"] + "tmp/"
    # config["dir_variants"] = config["dir_data"] + "variants/"
    # config["dir_logs"] = config["dir_data"] + "logs/"


    config["software_version_info"] = {}

    return False if error else config


config["teches"] = ["ILM", "BGI", "MGI", "CCS", "HiFi", "CLR", "ONTUL", "ONT"]
config["tech_short_reads"] = ["ILM", "BGI", "MGI"]
config["tech_long_reads"] = ["CCS", "HiFi", "CLR", "ULONT", "ONT"]
config["required_software"] = ["samtools"]
new_config = check_config_file(config)
if new_config:
    config.update(new_config)
else:
    print(f"[Error] The configure file(s) for this project is incorrect. "
          f"Please make the necessary changes according to the prompts")
    exit(-1)


rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    threads: get_run_threads("samtools_index")
    run:
        samtools = config["software"]["samtools"]
        shell("{samtools} index -@ {threads} {input}")

rule tabix_vcf:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    threads: get_run_threads("__default__")
    run:
        tabix = config["software"]["tabix"]
        shell("{tabix} -@ {threads} {input}")


def get_chroms_raw_vcf_merged_vcfs_cohort(wildcards):
    chroms = config["refs"][wildcards.ref_name]["available_chrom"]
    vcfs = []
    vcfs_idx = []
    for chrom in chroms:
        if wildcards.caller in ["hipstr"] and chrom in ["chrM"]: continue
        vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/chroms/{wildcards.cohort}.{wildcards.ref_name}.{wildcards.tech}."
                                         f"{wildcards.caller}.{chrom}.{wildcards.suffix}.raw.vcf.gz")
        vcfs_idx.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/chroms/{wildcards.cohort}.{wildcards.ref_name}.{wildcards.tech}.{wildcards.caller}.{chrom}.{wildcards.suffix}.raw.vcf.gz.tbi")
    return {"vcf": vcfs, "vcf_idx": vcfs_idx}


rule cohort_vcf_chrom_concat:
    input:
        unpack(get_chroms_raw_vcf_merged_vcfs_cohort)
    output:
        config["dir_data"] + "variants_raw/{cohort}/{caller}/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.raw.vcf.gz"
    log:
        config["dir_data"] + "variants_raw/{cohort}/{caller}/logs/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.chrom_concat.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/{caller}/logs/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.chrom_concat.rtime.tsv"
    wildcard_constraints:
        caller="varscan|bcftools|GATK_HC|pindel|manta|hipstr|gangstr|pbsv|freebayes|strelka|arcsv|hipstr|ExpansionHunter"
    threads: get_run_threads("SNVIndel_chrom_concat")
    run:
        bcftools = config["software"]["bcftools"]
        shell("{bcftools} concat -a  -o {output} -Oz --threads {threads} {input.vcf}  2>{log} 1>{log} ")


# def get_chroms_raw_vcf_merged_vcfs_sample(wildcards):
#     chroms = config["refs"][wildcards.ref_name]["available_chrom"]
#     vcfs = []
#     vcfs_idx = []
#     for chrom in chroms:
#         if wildcards.caller in ["longtr", "gangstr"] and chrom in ["chrM"]: continue
#         vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/samples/chroms/"
#                                          f"{wildcards.cohort}.{wildcards.sample}.{wildcards.ref_name}.{wildcards.tech}.{wildcards.caller}.{chrom}.{wildcards.suffix}.raw.vcf.gz")
#         vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/samples/chroms/"
#                                          f"{wildcards.cohort}.{wildcards.sample}.{wildcards.ref_name}.{wildcards.tech}.{wildcards.caller}.{chrom}.{wildcards.suffix}.raw.vcf.gz.tbi")
#
#     # vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/chroms/{wildcards.cohort}.{wildcards.ref_name}.{wildcards.tech}.{wildcards.caller}.{chrom}.{wildcards.suffix}.raw.vcf.gz")
#     # vcfs_idx.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/chroms/{wildcards.cohort}.{wildcards.ref_name}.{wildcards.tech}.{wildcards.caller}.{chrom}.{wildcards.suffix}.raw.vcf.gz.tbi")
#     return {"vcf": vcfs, "vcf_idx": vcfs_idx}


# rule sample_vcf_chrom_concat:
#     input:
#         unpack(get_chroms_raw_vcf_merged_vcfs_sample)
#     output:
#         # config["dir_data"] + "variants_raw/{cohort}/{caller}/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.raw.vcf.gz"
#         config["dir_data"] + "variants_raw/{cohort}/{caller}/samples/{cohort}.{sample}.{ref_name}.{tech}.{caller}.{suffix}.raw.vcf.gz",
#     log:
#         config["dir_data"] + "variants_raw/{cohort}/{caller}/logs/{cohort}.{sample}.{ref_name}.{tech}.{caller}.{suffix}.chrom_concat.log"
#     benchmark:
#         config["dir_data"] + "variants_raw/{cohort}/{caller}/logs/{cohort}.{sample}.{ref_name}.{tech}.{caller}.{suffix}.chrom_concat.rtime.tsv"
#     wildcard_constraints:
#         caller="longtr"
#     threads: get_run_threads("SNVIndel_chrom_concat")
#     run:
#         bcftools = config["software"]["bcftools"]
#         shell("{bcftools} concat  -o {output} -Oz --threads {threads} {input.vcf}  2>{log} 1>{log} ")


def get_chroms_raw_vcf_merged_vcfs_sample_new(wildcards):
    chroms = config["refs"][wildcards.ref_name]["available_chrom"]
    vcfs = []
    vcfs_idx = []
    for chrom in chroms:
        if wildcards.caller in ["longtr", "gangstr","ExpansionHunter"] and chrom in ["chrM"]: continue
        # config["dir_data"] + "variants_raw/{cohort}/longshot/samples/{cohort}.{sample}.{ref_name}.{tech}/{chrom}.longshot.SNVIndel.raw.vcf.gz",
        vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/samples/"
                                         f"{wildcards.cohort}.{wildcards.sample}.{wildcards.ref_name}.{wildcards.tech}/{chrom}.{wildcards.caller}.{wildcards.suffix}.raw.vcf.gz")
        vcfs_idx.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/samples/"
                                             f"{wildcards.cohort}.{wildcards.sample}.{wildcards.ref_name}.{wildcards.tech}/{chrom}.{wildcards.caller}.{wildcards.suffix}.raw.vcf.gz.tbi")

    # vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/chroms/{wildcards.cohort}.{wildcards.ref_name}.{wildcards.tech}.{wildcards.caller}.{chrom}.{wildcards.suffix}.raw.vcf.gz")
    # vcfs_idx.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/chroms/{wildcards.cohort}.{wildcards.ref_name}.{wildcards.tech}.{wildcards.caller}.{chrom}.{wildcards.suffix}.raw.vcf.gz.tbi")
    return {"vcf": vcfs, "vcf_idx": vcfs_idx}


rule sample_vcf_chrom_concat_new:
    input:
        unpack(get_chroms_raw_vcf_merged_vcfs_sample_new)
    output:
        # config["dir_data"] + "variants_raw/{cohort}/{caller}/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.raw.vcf.gz"
        config["dir_data"] + "variants_raw/{cohort}/{caller}/samples/{cohort}.{sample}.{ref_name}.{tech}.{caller}.{suffix}.raw.vcf.gz",
    log:
        config["dir_data"] + "variants_raw/{cohort}/{caller}/logs/{cohort}.{sample}.{ref_name}.{tech}.{caller}.{suffix}.chrom_concat.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/{caller}/logs/{cohort}.{sample}.{ref_name}.{tech}.{caller}.{suffix}.chrom_concat.rtime.tsv"
    wildcard_constraints:
        caller="longshot|clair3|nanocaller|nanovar|longtr|hipstr|gangstr|ExpansionHunter"
    threads: get_run_threads("SNVIndel_chrom_concat")
    run:
        bcftools = config["software"]["bcftools"]
        shell("{bcftools} concat  -o {output} -Oz --threads {threads} {input.vcf}  2>{log} 1>{log} ")

def get_merge_samples_input(wildcards):
    vcfs = []
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        # bams.append(path_bam)
        # bais.append(path_bam + ".bai")
        vcfs.append(config["dir_data"] + f"variants_raw/{wildcards.cohort}/{wildcards.caller}/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.{wildcards.caller}.{wildcards.suffix}.vcf.gz")
    return vcfs

rule merge_samples_vcfs:
    input:
        unpack(get_merge_samples_input)
    output:
        config["dir_data"] + "variants_raw/{cohort}/{caller}/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.vcf.gz"
    wildcard_constraints:
        caller="longtr|hipstr|gangstr"
    run:
        shell("")
