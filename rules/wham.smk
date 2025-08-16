rule whamg_call:
    input:
        unpack(get_samples_bam),
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        vcf=config["dir_data"] + "variants_raw/{cohort}/whamg/{cohort}.{ref_name}.{tech}.whamg.SV.raw.vcf.gz"
    params:
        extra="",
    log:
        config["dir_data"] + "variants_raw/{cohort}/whamg/logs/{cohort}.{ref_name}.{tech}.whamg.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/whamg/logs/{cohort}.{ref_name}.{tech}.whamg.rtime.tsv"

    threads: get_run_threads("whamg_call")
    run:
        whamg = config["software"]["whamg"]
        bcftools = config["software"]["bcftools"]
        exclude_contigs = []
        chroms = config["refs"][wildcards.ref_name]["available_chrom"]
        for line in open(f"{input.ref}" + ".fai"):
            chrom = line.split("\t")[0]
            if chrom not in chroms:
                exclude_contigs.append(chrom)
        exclude_contig_str = ",".join([f"{i}" for i in exclude_contigs])
        bam_str = ",".join([i for i in input.bams])
        out_vcf = f"{output.vcf}"[:-3]
        shell("export EXCLUDE={exclude_contig_str} && "
              "{whamg} -x {threads} -e $EXCLUDE -a {input.ref} -f {bam_str} > {out_vcf}  2> {log} ")
        shell("{bcftools} view -o {output.vcf} -Oz  {out_vcf}")

rule wham_call:
    input:
        unpack(get_samples_bam),
        exclude=lambda wildcards: config["refs"][wildcards.ref_name]["exclude"],
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
    output:
        vcf=config["dir_data"] + "variants_raw/{cohort}/wham/{cohort}.{ref_name}.{tech}.wham.SV.raw.vcf.gz"
    params:
        extra="",
    log:
        config["dir_data"] + "variants_raw/{cohort}/wham/logs/{cohort}.{ref_name}.{tech}.wham.log"
    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/wham/logs/{cohort}.{ref_name}.{tech}.wham.rtime.tsv"

    threads: get_run_threads("wham_call")
    run:
        wham = config["software"]["wham"]
        bcftools = config["software"]["bcftools"]
        bam_str = ",".join([i for i in input.bams])
        out_vcf = f"{output.vcf}"[:-3]
        shell("{wham} -x {threads} -e {input.exclude} -f {input.ref} -t {bam_str} > {out_vcf}  2> {log} ")
        shell("{bcftools} view -o {output.vcf} -Oz  {out_vcf}")
