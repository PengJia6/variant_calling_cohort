#TODO process exclude cnvnator
rule delly_call:
    input:
        unpack(get_samples_bam),
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        vcf=config["dir_data"] + "variants_raw/{cohort}/delly/{cohort}.{ref_name}.{tech}.delly.SV.raw.vcf.gz",
        bcf=config["dir_data"] + "variants_raw/{cohort}/delly/{cohort}.{ref_name}.{tech}.delly.SV.raw.bcf",
        bcf_gt=config["dir_data"] + "variants_raw/{cohort}/delly/{cohort}.{ref_name}.{tech}.delly.SV.4gt.bcf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/delly/logs/{cohort}.{ref_name}.{tech}.delly.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/delly/logs/{cohort}.{ref_name}.{tech}..delly.rtime.tsv"
    threads: get_run_threads("__default_single_threads")
    run:
        delly = config["software"]["delly"]
        bcftools = config["software"]["bcftools"]
        if f"{wildcards.tech}".upper() in ["ONT", "ONTUL", "ULONT", "HIFI", "CCS", "CLR"]:
            shell("{delly} lr -g {input.ref}  -o {output.bcf_gt} {input.bams} 2>{log} ")
            shell("{delly} lr -g {input.ref}  -v {output.bcf_gt} -o {output.bcf}  {input.bams} 2>>{log} ")
            shell("{bcftools} sort -Oz  -o {output.vcf} {output.bcf}")

        elif f"{wildcards.tech}".upper() in ["NGS", "ILM", "MGI", "BGI"]:
            # shell("{delly} call -g {input.ref} -x {input.excl} -o {output.bcf} {input.bam} 2>{log} 1>{log}")
            shell("{delly} call -g {input.ref}  -o {output.bcf_gt} {input.bams} 2>{log} ")
            shell("{delly} call -g {input.ref}  -v {output.bcf_gt} -o {output.bcf}  {input.bams} 2>>{log} ")
            shell("{bcftools} sort -Oz  -o {output.vcf} {output.bcf}")
        else:
            print("NO available model provided")
            exit()
