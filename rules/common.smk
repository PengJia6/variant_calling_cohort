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


def get_one_sample_fq_R1(wildcards):
    return config["samples"][wildcards.cohort][wildcards.sample]["path"][wildcards.sample_type][wildcards.tech][
        wildcards.subsample]["R1"]


def get_one_sample_fq_R2(wildcards):
    return config["samples"][wildcards.cohort][wildcards.sample]["path"][wildcards.sample_type][wildcards.tech][
        wildcards.subsample]["R2"]


def get_one_sample_fq(wildcards):
    return config["samples"][wildcards.cohort][wildcards.sample]["path"][wildcards.sample_type][wildcards.tech][
        wildcards.subsample]


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return f"'@RG\\tID:{wildcards.subsample}\\tLB:{wildcards.subsample}\\tSM:{wildcards.sample}_{wildcards.sample_type}\\tPL:{wildcards.tech}' "


def check_samples(conf_samples):
    error = False
    file_samples = open(conf_samples,"r")
    samples_info = yaml.safe_load(file_samples)
    config["samples"] = samples_info["samples"]
    file_samples.close()
    # if "samples" not in config:
    # samples_info = config["samples"]
    # print(samples_info)
    cohort_names = []
    sample_names = []
    sample_types = []
    tech_names = []
    lib_names = []

    for cohort, cohort_info in samples_info["samples"].items():
        cohort_names.append(cohort)
        for sample_name, info_t in cohort_info.items():
            sample_names.append(sample_name)
            if "path" not in info_t:
                print("[Error] There is not 'path' info in sample configure file!")
                error = True

            for sample_type, info_tt in info_t["path"].items():
                sample_types.append(sample_type)
                if len(info_tt) < 1:
                    print(f"[Error] There is not available reads in {sample_type} sample of {sample_name} !")
                    error = True

                for tech, info in info_tt.items():
                    tech_names.append(tech)
                    if len(info) < 1:
                        print(f"[Error] There is not available reads in sample {tech} tech of {sample_type} sample of {sample_name}")
                        error = True
                    if ("paras" not in config) or (tech not in config["paras"]) or "aligners" not in config["paras"][tech]:
                        print(f"[Error] error setting of alignment paras for {tech} reads alignment in configure file. e.g.\n"
                              f"-----------------------\n...\n"
                              f"paras:\n\t{tech}:\n\t\taligners:\n\t\t\t- xxx\n\t\t\t- yyy、n"
                              f"-----------------------\n")
                        error = True
                    else:
                        for aligner in config["paras"][tech]["aligners"]:
                            config["required_software"].append(aligner)
                        if "mark_duplication" in config["paras"][tech] and config["paras"][tech]["mark_duplication"]:
                            config["required_software"].append("bammarkduplicates2")
                        if "bam_left_align" in config["paras"][tech] and config["paras"][tech]["bam_left_align"]:
                            config["required_software"].append("gatk")
                    if tech in config["tech_short_reads"]:
                        for sub_name, one_path in info.items():
                            # print(sub_name, one_path)
                            if "R1" not in one_path or "R2" not in one_path:
                                print(
                                    f"[Error] Please provide paired-end path of fastq file for {sub_name} of {sample_name}")
                            if not (one_path["R1"].endswith("fastq.gz") or one_path["R1"].endswith("fq.gz") or
                                    one_path["R1"].endswith("fastq") or one_path["R1"].endswith("fq") or
                                    one_path["R1"].endswith("fasta") or one_path["R1"].endswith("fa") or
                                    one_path["R1"].endswith("fasta.gz") or one_path["R1"].endswith("fa.gz")
                            ):
                                print("[Error] The path {} not end with fastq or fasta!".format(
                                    one_path["R1"]))
                                error = True
                            if not (one_path["R2"].endswith("fastq.gz") or (one_path["R2"].endswith("fq.gz")) or
                                    one_path["R2"].endswith("fastq") or (one_path["R2"].endswith("fq")) or
                                    one_path["R2"].endswith("fasta.gz") or (one_path["R2"].endswith("fa.gz")) or
                                    one_path["R2"].endswith("fasta") or (one_path["R2"].endswith("fa"))
                            ):
                                print("[Error] The path {} not end with fastq.gz, fq.gz, bam, or cram !".format(
                                    one_path["R2"]))
                                error = True
                            if one_path["R1"] == one_path["R2"]:
                                print(f"[Error] The R1 path ({one_path['R1']}) same as R2 path ({one_path['R2']})!!!!")
                                error = True

                    elif tech in config["tech_long_reads"]:
                        # this_sub_samples = {}
                        for lib, one_path in info.items():
                            lib_names.append(lib)
                            if not (one_path.endswith("fastq.gz") or one_path.endswith("fq.gz") or
                                    one_path.endswith("fa") or one_path.endswith("fasta") or
                                    one_path.endswith("fa.gz") or one_path.endswith("fasta.gz") or
                                    one_path.endswith("fastq") or one_path.endswith("fq")):
                                print(f"[Error] The path '{one_path}' not end with fastq.gz, fq.gz, fastq, fq, fa, fa.gz, fasta, or fasta.gz !")
                                error = True

    config["name_space"] = {"cohort_names": list(set(cohort_names)),
                            "sample_names": list(set(sample_names)),
                            "sample_types": list(set(sample_types)),
                            "tech_names": list(set(tech_names)),
                            "lib_names": list(set(lib_names))}

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
        config["path_samples_config"] = "confs/alignment_samples.yaml"
    if not os.path.exists(config["path_samples_config"]):
        print(
            f"[Error] The config file {config['path_samples_config']} for samples ('path_samples_config' in config.yaml) is not available! "
            "Please provide a config file (yaml format) for the reads of the samples in this project! "
            "You can setting the config file according the template xxxx")  # TODO add the link
        error = True
    else:
        check_samples(config["path_samples_config"])
        # print(samples_info)

        if not check_samples(config["path_samples_config"]):
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
