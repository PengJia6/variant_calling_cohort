#TODO process exclude cnvnator
rule delly_call:
    input:
        bam=config["dir_aligned_reads"] + "{prefix}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{prefix}.{ref_name}.{suffix}.bam.bai",
        ref=config["dir_ref"] + "{ref_name}.fasta",
        excl="/data/DATA/Reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38.d1.vd1.exclude.cnvnator.bed"
    output:
        bcf=config["dir_variants"] + "{prefix}.{ref_name}.{suffix}.delly.discover.bcf",
    log:
        config["dir_variants"] + "{prefix}.{ref_name}.{suffix}..delly.discover.log",
    benchmark:
        config["dir_variants"] + "{prefix}.{ref_name}.{suffix}..delly.discover.rtime.tsv"

    threads: get_run_threads("delly_gt")
    run:
        delly = config["software"]["delly"]
        shell("{delly} call -g {input.ref} -x {input.excl} -o {output.bcf} {input.bam} 2>{log} 1>{log}")


rule delly_gt:
    input:
        bam=config["dir_aligned_reads"] + "{prefix}.{ref_name}.{suffix}.bam",
        bai=config["dir_aligned_reads"] + "{prefix}.{ref_name}.{suffix}.bam.bai",
        ref=config["dir_ref"] + "{ref_name}.fasta",
        excl="/data/DATA/Reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38.d1.vd1.exclude.cnvnator.bed",
        bcf=config["dir_variants"] + "{prefix}.{ref_name}.{suffix}.delly.discover.bcf",
    output:
        bcf=config["dir_variants"] + "{prefix}.{ref_name}.{suffix}.delly.raw.bcf",
        vcf=config["dir_variants"] + "{prefix}.{ref_name}.{suffix}.delly.raw.vcf.gz",
    log:
        config["dir_variants"] + "{prefix}.{ref_name}.{suffix}.delly_gt.logs"
    benchmark:
        config["dir_variants"] + "{prefix}.{ref_name}.{suffix}.delly_gt.logs"
    threads: get_run_threads("delly_gt")
    run:
        bcftools = config["software"]["bcftools"]
        delly = config["software"]["delly"]
        shell("{delly} call -g {input.ref} -x {input.excl} -o {output.bcf} -v {input.bcf} {input.bam} 2>{log} 1>{log}")
        shell("{bcftools} view -Oz -o {output.vcf} {output.bcf}")
