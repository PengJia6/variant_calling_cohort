# ======================================================================================================================
# Project: HiFi_analysis_pipeline
# Script : common.py.smk
# Author : Peng Jia
# Date   :  2022/9/29
# Email  : pengjia@stu.xjtu.edu.cn
# Description: common functions for this pipeline
# ======================================================================================================================

import os.path

import yaml


def get_run_threads(rule_name):
    if rule_name in config["threads"]:
        return config["threads"][rule_name]["cpus"]
    else:
        return config["threads"]["__default__"]["cpus"]


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
        if "path" not in cohort_info or "info" not in cohort_info or "tech" not in cohort_info:
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
            config["refs"] = yaml.safe_load(open(config["path_ref_config"]))
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

    if "path_cluster_config" not in config:
        config["path_cluster_config"] = "confs/alignment_cluster.yaml"
    if not os.path.exists(config["path_cluster_config"]):
        print(f"[Error] The config file {config['path_cluster_config']} for resources"
              f" ('path_cluster_config' in config.yaml) is not available! "
              "Please provide a config file (yaml format) to ensure that each rule utilizes a sensible threads! "
              "You can setting the config file according the template xxxx")  # TODO add the link
        error = True
    else:
        config["threads"] = yaml.safe_load(open(config["path_cluster_config"]))

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

    config["dir_data"] = os.path.abspath("align_reads") if "dir_data" not in config else config["dir_data"]
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
    threads: get_run_threads("bwa_align")
    run:
        samtools = config["software"]["samtools"]
        shell("{samtools} index -@ {threads} {input}")

rule tabix_vcf:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    threads: get_run_threads("tabix_vcf")
    run:
        tabix = config["software"]["tabix"]
        shell("{tabix} -@ {threads} {input}")


def get_chroms_raw_vcf_merged_vcfs(wildcards):
    chroms = config["refs"][wildcards.ref_name]["available_chrom"]

    return expand(config["dir_data"] + "{cohort}/{caller}/chroms/{cohort}.{ref_name}.{tech}.{caller}.{chrom}.{suffix}.raw.vcf.gz",
        cohort=wildcards.cohort,caller=wildcards.caller,ref_name=wildcards.ref_name,tech=wildcards.tech,chrom=chroms,suffix=wildcards.suffix)


rule SNVIndel_chrom_concat:
    input:
        get_chroms_raw_vcf_merged_vcfs
    output:
        config["dir_data"] + "{cohort}/{caller}/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.raw.vcf.gz"
    log:
        config["dir_data"] + "{cohort}/logs/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.SNVIndel_chrom_concat.log"
    benchmark:
        config["dir_data"] + "{cohort}/logs/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.SNVIndel_chrom_concat.rtime.tsv"
    wildcard_constraints:
        caller="varscan|bcftools|GATK_HC"
    threads: get_run_threads("varscan_call_snp_indel")
    run:
        bcftools = config["software"]["bcftools"]
        shell("{bcftools} concat  -o {output} -Oz --threads {threads} {input}  2>{log} 1>{log} ")
