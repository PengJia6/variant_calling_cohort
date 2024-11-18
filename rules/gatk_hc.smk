# ======================================================================================================================
# Project: variant_calling_cohort
# Script : varscan.smk TODO check 
# Author : Peng Jia
# Date   :  2024/11/18
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: TODO
# ======================================================================================================================

rule gatk_hc_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    # ref_dict= lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        gvcf=config["dir_data"] + "{cohort}/GATK_HC/samples/{sample}_{tech}/{cohort}.{sample}.{ref_name}.{tech}.GATK_HC.{chrom}.SNVIndel.g.vcf.gz"
    # =config["dir_variants"] + "gatk/gatk_details/{sample}/{sample}.{prefix}.{contig}.gvcf.gz"
    log:
        config["dir_data"] + "{cohort}/GATK_HC/samples/{sample}_{tech}_logs/{cohort}.{sample}.{ref_name}.{tech}.GATK_HC.{chrom}.SNVIndel.gatk_hc.log"
    threads: get_run_threads("gatk_hc_call")
    benchmark:
        config["dir_data"] + "{cohort}/GATK_HC/samples/{sample}_{tech}_logs/{cohort}.{sample}.{ref_name}.{tech}.GATK_HC.{chrom}.SNVIndel.gatk_hc.rtime.tsv"

    params:
        extra="",
        java_options="",
        regions="",
        dbsnp=[],
    run:
        gatk = config["software"]["gatk"]
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{gatk} {java_opt} HaplotypeCaller {params.extra} --minimum-mapping-quality 8  "
              " -R {input.ref} -ERC GVCF -L {wildcards.chrom} -I {input.bam} -O {output.gvcf}"
              " 2>{log} 1>{log}")


#

def get_conbine_gvcf_input(wildcards):
    vcfs = []
    # for sample_name, info_t in config["samples_info"][f"{wildcards.cohort}"].items():
    #     for sample_type, info_tt in info_t["path"].items():
    #         # print(wildcards.cohort,sample_type,sample_name,wildcards.this_prefix)
    #         # for tech, info in info_tt.items():
    #         vcfs.append(config["dir_variants"] + f"{wildcards.cohort}/{sample_name}/{sample_type}/{wildcards.cohort}.{sample_name}.{sample_type}.{wildcards.ref_name}.{wildcards.this_prefix}.gatk.{wildcards.contig}.raw.g.vcf.gz")
    bams = []
    bais = []
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        bams.append(path_bam)
        bais.append(path_bam + ".bai")
        vcfs.append(config["dir_data"] + f"{wildcards.cohort}/GATK_HC/samples/{sample}_{wildcards.tech}/{wildcards.cohort}.{sample}.{wildcards.ref_name}.{wildcards.tech}.GATK_HC.{wildcards.chrom}.SNVIndel.g.vcf.gz",)
    return {"bams": bams, "bais": bais,
            "ref": config["refs"][wildcards.ref_name]["fasta"],
            "ref_idx": config["refs"][wildcards.ref_name]["fasta"] + ".fai",
            "gvcfs": vcfs
            }
    return vcfs


#

rule gatk_combine_gvcf:
    input:
        unpack(get_conbine_gvcf_input)
    output:
        vcf=config["dir_data"] + "{cohort}/GATK_HC/chroms/{cohort}.{ref_name}.{tech}.GATK_HC.{chrom}.SNVIndel.g.vcf.gz",
    params:
        extra="",
        java_options=""
    threads: get_run_threads("gatk_combine_gvcf")
    log:
        config["dir_data"] + "{cohort}/GATK_HC/logs/{cohort}.{ref_name}.{tech}.GATK_HC.{chrom}.SNVIndel.gatk_combine_gvcf.log"
    benchmark:
        config["dir_data"] + "{cohort}/GATK_HC/logs/{cohort}.{ref_name}.{tech}.GATK_HC.{chrom}.SNVIndel.gatk_combine_gvcf.rtime.tsv"

    run:
        gatk = config["software"]["gatk"]
        inputs = " ".join([("-V " + f) for f in input.gvcfs])
        input_bams = " ".join(["-I {i}" for i in input.bams])

        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{gatk} {java_opt} CombineGVCFs {params.extra} "
              " {inputs} {input_bams} -R {input.ref} -O {output} 2>{log} 1>{log}")


def get_gatk_genotype_input(wildcards):
    bams = []
    bais = []
    for sample, path_bam in config["samples"][wildcards.cohort]["path"][wildcards.tech].items():
        bams.append(path_bam)
        bais.append(path_bam + ".bai")
    return {"bams": bams, "bais": bais,
            "ref": config["refs"][wildcards.ref_name]["fasta"],
            "ref_idx": config["refs"][wildcards.ref_name]["fasta"] + ".fai",
            "vcf": config["dir_data"] + f"{wildcards.cohort}/GATK_HC/chroms/{wildcards.cohort}.{wildcards.ref_name}.{wildcards.tech}.GATK_HC.{wildcards.chrom}.SNVIndel.g.vcf.gz",
            }


rule gatk_genotype:
    input:
        unpack(get_gatk_genotype_input)
    # vcf=config["dir_variants"] + "{cohort}.{ref_name}.{this_prefix}.gatk.{contig}.g.vcf.gz"

    output:
        vcf=config["dir_data"] + "{cohort}/GATK_HC/chroms/{cohort}.{ref_name}.{tech}.GATK_HC.{chrom}.SNVIndel.raw.vcf.gz",
    # vcf=config["dir_variants"] + "gatk/gatk_details/contigs/" + config["project"] + ".{prefix}.{contig}.vcf.gz"
    params:
        extra="",
        java_options=""
    threads: get_run_threads("gatk_genotype")
    log:
        config["dir_data"] + "{cohort}/GATK_HC/logs/{cohort}.{ref_name}.{tech}.varscan.{chrom}.SNVIndel.GATK_HC.log"
    benchmark:
        config["dir_data"] + "{cohort}/GATK_HC/logs/{cohort}.{ref_name}.{tech}.varscan.{chrom}.SNVIndel.GATK_HC.rtime.tsv"
    run:
        gatk = config["software"]["gatk"]
        input_bams = " ".join(["-I {i}" for i in input.bams])
        java_opt = "" if len(params.java_options) < 2 else " --java-options {params.java_options} "
        shell("{gatk} {java_opt}  GenotypeGVCFs {params.extra} "
              " -R {input.ref} -V {input.vcf} {input_bams}  -O {output.vcf} 2>{log} 1>{log}")
