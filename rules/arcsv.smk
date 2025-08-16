rule arcsv_call:
    input:
        unpack(get_samples_bam),
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        gap= lambda wildcards: config["refs"][wildcards.ref_name]["gap"],

    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/arcsv/chroms/{cohort}.{ref_name}.{tech}.arcsv.{chrom}.SV.raw.vcf.gz"
    # prefix=config["dir_data"] + "variants_raw/{cohort}/gridss/chroms/{cohort}.{ref_name}.{tech}.{chrom}.gridss/runWorkflow.py"
    params:
        extra="",
    # dp=5
    log:
        config["dir_data"] + "variants_raw/{cohort}/arcsv/logs/{cohort}.{ref_name}.{tech}.{chrom}.arcsv.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/arcsv/logs/{cohort}.{ref_name}.{tech}.{chrom}.arcsv.rtime.tsv"
    threads: get_run_threads("arcsv_call")

    run:
        arcsv = config["software"]["arcsv"]
        bcftools = config["software"]["bcftools"]
        prefix = str(output.vcfgz)[:-7]
        shell("mkdir -p {prefix}")
        bams_str = ",".join([f"{i}" for i in input.bams])
        vcf_path = f"{output.vcfgz}"[:-6]+"unsort.vcf.gz"
        # vcf_path_sort=f"{output.vcfgz}"[:-3]
        shell("{arcsv} call  -i {bams_str} -r {wildcards.chrom} -R {input.ref}  -G {input.gap} --overwrite "
              "-o {prefix} 2>{log} 1>{log}")
        #
        # shell("{bcftools} view -Oz -o {vcf_path} {prefix}/arcsv_out.vcf")

        shell("{bcftools} sort  -o {output.vcfgz} -Oz {prefix}/arcsv_out.vcf")

# shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")

# /data/home/lbjia/miniconda3/envs/default/bin/bcftools
