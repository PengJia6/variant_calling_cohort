# ======================================================================================================================
# Project: variant_calling_cohort
# Script : varscan.smk TODO check 
# Author : Peng Jia
# Date   :  2024/11/18
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: TODO
# ======================================================================================================================
rule deepvariant_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcf_gz=config["dir_data"] + "{cohort}/deepvariant/samples/{cohort}.{sample}.{ref_name}.{tech}.deepvariant.SNVIndel.vcf.gz",
        gvcf_gz=config["dir_data"] + "{cohort}/deepvariant/samples/{cohort}.{sample}.{ref_name}.{tech}.deepvariant.SNVIndel.g.vcf.gz"
    # gvcf_gz=config["dir_variants"] + "dv/dv_details/{sample}/{sample}.{prefix}.dv.raw.g.vcf.gz"
    log:
        config["dir_data"] + "{cohort}/deepvariant/samples/{cohort}.{sample}.{ref_name}.{tech}.deepvariant.SNVIndel.log",

    benchmark:
        config["dir_data"] + "{cohort}/deepvariant/samples/{cohort}.{sample}.{ref_name}.{tech}.deepvariant.SNVIndel.rtime.tsv",
    threads: get_run_threads("deepvariant_call")
    run:
        dir_tmp = str(output.vcf_gz).rstrip(".vcf.gz") + "_tmp"
        file_tmp = dir_tmp.split("/")[-1]
        shell("mkdir -p " + dir_tmp)
        bam_dir = "/".join(str(input.bam).split("/")[:-1])
        bam_file = str(input.bam).split("/")[-1]
        ref_dir = "/".join(str(input.ref).split("/")[:-1])
        ref_file = str(input.ref).split("/")[-1]
        output_dir = "/".join(str(output.vcf_gz).split("/")[:-1])
        output_file = str(output.vcf_gz).split("/")[-1].rstrip(".vcf.gz")
        bcftools = config["software"]["bcftools"]
        if wildcards.tech in ["ILM", "BGI", "WGS", "NGS"]:
            this_tech = "WGS"
        elif wildcards.tech in ["CCS", "HiFi", "PACBIO"]:
            this_tech = "PACBIO"
        elif wildcards.tech in ["WES"]:
            this_tech = "WES"
        else:
            exit()
        shell('docker run '
              '-v "{bam_dir}":"/input" '
              '-v "{ref_dir}":"/ref" '
              '-v "{output_dir}":"/output" '
              'google/deepvariant:1.6.0 /opt/deepvariant/bin/run_deepvariant '
              '--model_type={this_tech} '
              '--ref=/ref/{ref_file} '
              '--reads=/input/{bam_file} '
              '--output_vcf=/output/{output_file}.vcf '
              '--output_gvcf=/output/{output_file}.g.vcf '
              '--num_shards={threads} '
              '--make_examples_extra_args min_mapping_quality=1,keep_supplementary_alignments=true '
              '--intermediate_results_dir /output/{file_tmp} 1>>{log} 2>>{log}')
        shell("{bcftools} view -Oz -o {output.vcf_gz} {output_dir}/{output_file}.vcf 2>>{log} 1>>{log}")
        shell("{bcftools} view -Oz -o {output.gvcf_gz} {output_dir}/{output_file}.g.vcf 2>>{log} 1>>{log}")


def get_vcfs_for_deepvariant_merge_vcf(wildcards):
    gvcfs = []
    vcfs = []
    # for cohort, cohort_info in.items():
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        # for sample_name, info_t in config["samples_info"][f"{wildcards.cohort}"].items():
        #     for sample_type, info_tt in info_t["path"].items():
        # print(wildcards.cohort,sample_type,sample_name,wildcards.this_prefix)
        # for tech, info in info_tt.items():
        vcfs.append(config["dir_data"] + f"{wildcards.cohort}/deepvariant/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.deepvariant.SNVIndel.vcf.gz")
        gvcfs.append(config["dir_data"] + f"{wildcards.cohort}/deepvariant/samples/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.deepvariant.SNVIndel.g.vcf.gz")
    return {"gvcf_gz": gvcfs, "vcf_gz": vcfs}


#

rule deepvariant_merge_vcf:
    input:
        unpack(get_vcfs_for_deepvariant_merge_vcf)
    output:
        vcf=config["dir_data"] + "{cohort}/deepvariant/chroms/{cohort}.{ref_name}.{tech}.deepvariant.{chrom}.SNVIndel.raw.vcf.gz",
    params: ""
    log:
        config["dir_data"] + "{cohort}/deepvariant/logs/{cohort}.{ref_name}.{tech}.{chrom}.deepvariant.merge.log"
    benchmark:
        config["dir_data"] + "{cohort}/deepvariant/logs/{cohort}.{ref_name}.{tech}.{chrom}.deepvariant.merge.rtime.tsv"
    run:
        if len(input.gvcf_gz) == 1:
            shell("cp {input.vcf_gz} {output.vcf}")
        else:
            input_vcfs = []
            input_dir = "/".join(input.gvcf_gz[0].split("/")[:-4])
            for item in input.gvcf_gz:
                input_vcfs.append("/".join(item.split("/")[-4:]))
            inputs_str = " ".join(["/input/" + item for item in input_vcfs])
            bcftools = config["software"]["bcftools"]

            shell("docker run -v {input_dir}:/input quay.io/mlin/glnexus:v1.2.7 /usr/local/bin/glnexus_cli "
                  "--config DeepVariantWGS "
                  "{inputs_str} |"
                  "{bcftools} view -Oz -o {output.vcf}")
