rule gridss_call:
    input:
        unpack(get_samples_bam),
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/gridss/{cohort}.{ref_name}.{tech}.gridss.SV.raw.vcf.gz"
    # prefix=config["dir_data"] + "variants_raw/{cohort}/gridss/chroms/{cohort}.{ref_name}.{tech}.{chrom}.gridss/runWorkflow.py"
    params:
        extra="",
    # dp=5
    log:
        config["dir_data"] + "variants_raw/{cohort}/gridss/logs/{cohort}.{ref_name}.{tech}.gridss.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/gridss/logs/{cohort}.{ref_name}.{tech}.gridss.rtime.tsv"
    threads: get_run_threads("gridss_call")
    # config["dir_data"] + "variants_raw/{cohort}/{caller}/{cohort}.{ref_name}.{tech}.{caller}.{suffix}.raw.vcf.gz"
    run:
        gridss = config["software"]["gridss"]
        bcftools = config["software"]["bcftools"]
        prefix = str(output.vcfgz)[:-7]
        shell("mkdir -p {prefix}")
        bams_str = " ".join([f"{i}" for i in input.bams])
        vcf_path = f"{output.vcfgz}"[:-3]
        # vcf_path_sort=f"{output.vcfgz}"[:-3]
        gridss_preifx = "/".join(f"{gridss}".split("/")[:-1])
        shell("cd {prefix} && "
              "export PATH={gridss_preifx}:$PATH && "
              " {gridss}  -r {input.ref} -t {threads} -w . {bams_str} -o {vcf_path} 2>{log} 1>{log}")
        shell("{bcftools} sort  -o {output.vcfgz} -Oz {vcf_path}")
# shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")

# /data/home/lbjia/miniconda3/envs/default/bin/bcftools
