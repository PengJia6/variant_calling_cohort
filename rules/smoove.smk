# Note: the duphold should be update to the latest version
rule smoove_call:
    input:
        unpack(get_samples_bam),
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        vcf=config["dir_data"] + "variants_raw/{cohort}/smoove/{cohort}.{ref_name}.{tech}.smoove.SV.raw.vcf.gz"
    # prefix=config["dir_data"] + "variants_raw/{cohort}/manta/chroms/{cohort}.{ref_name}.{tech}.{chrom}.manta/runWorkflow.py"
    params:
        extra="",
    log:
        config["dir_data"] + "variants_raw/{cohort}/smoove/logs/{cohort}.{ref_name}.{tech}.smoove.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/smoove/logs/{cohort}.{ref_name}.{tech}..smoove.rtime.tsv"

    threads: get_run_threads("smoove_call")
    run:
        smoove = config["software"]["smoove"]
        bcftools = config["software"]["bcftools"]
        prefix = str(output.vcf)[:-7]
        lumpy_prefix = "/".join(str(smoove).split("/")[:-1])
        import os

        if os.path.exists(f"{prefix}"):
            shell("rm -rf {prefix}")
        shell("mkdir -p {prefix}")
        shell("export PATH={lumpy_prefix}:$PATH && "
              "{smoove} call -x -d --genotype --name {wildcards.cohort} --outdir {prefix} --exclude {input.exclude} "
              "--fasta {input.ref} -p {threads} {input.bams} 2>{log} 1>{log}")
        shell("cp {prefix}/{wildcards.cohort}-smoove.genotyped.vcf.gz {output.vcf}")
        shell("touch {output.vcf}")
