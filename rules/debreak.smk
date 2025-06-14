# ======================================================================================================================
# Project: variant_calling_cohort
# Script : debreak.smk
# Author : Peng Jia
# Date   :  2025/04/17
# Email  : pengjia@xjtu.edu.cn; pengjia1110@163.com
# Description: debreak call
# ======================================================================================================================
rule debreak_call:
    input:
        bam=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample],
        bai=lambda wildcards: config["samples"][wildcards.cohort]["path"][wildcards.tech][wildcards.sample] + ".bai",
        ref=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"],
        ref_fai=lambda wildcards: config["refs"][wildcards.ref_name]["fasta"] + ".fai",
    # bam=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam",
    # bai=config["dir_aligned_reads"] + "{prefix}{ref_name}{suffix}.bam.bai",
    # ref=config["dir_ref"] + "{ref_name}.fasta",
    output:
        vcfgz=config["dir_data"] + "variants_raw/{cohort}/debreak/samples/{cohort}.{sample}.{ref_name}.{tech}.debreak.SV.raw.vcf.gz",
        vcf=config["dir_data"] + "variants_raw/{cohort}/debreak/samples/{cohort}.{sample}.{ref_name}.{tech}.debreak.SV.raw.vcf",
    # snf=config["dir_data"] + "{cohort}/sniffles/samples/{cohort}.{sample}.{ref_name}.{tech}.sniffles.SV.raw.snf",
    log:
        config["dir_data"] + "variants_raw/{cohort}/debreak/logs/{cohort}.{sample}.{ref_name}.{tech}.debreak.SV.log",

    benchmark:
        config["dir_data"] + "variants_raw/{cohort}/debreak/logs/{cohort}.{sample}.{ref_name}.{tech}.debreak.SV.rtime.tsv",
    threads: get_run_threads("debreak_call")
    run:
        workdir = str(output.vcfgz).rstrip(".vcf.gz") + "_tmp"
        debreak = config["software"]["debreak"]
        bcftools = config["software"]["bcftools"]
        minimap2 = config["software"]["minimap2"]
        wtdbg2 = config["software"]["wtdbg2"]
        prefix = ""
        for i in [debreak, minimap2, wtdbg2]:
            prefix += "/".join(i.split("/")[:-1])
        if os.path.exists(f"{workdir}"):
            shell("rm -rf {workdir}")

        shell("mkdir -p {workdir}")
        shell("date > {log}")
        shell("""
             echo  {debreak} -t {threads} --bam {input.bam} -o {workdir}/ --rescue_large_ins --rescue_dup --poa --ref {input.ref} --min_support 3 2>>{log} 1>>{log}  
                """)
        shell("""
        export PATH={prefix}$PATH && 
        {debreak} -t {threads} --bam {input.bam} -o {workdir}/ --rescue_large_ins --rescue_dup --poa --ref {input.ref} --min_support 3 2>>{log} 1>>{log}  
        """)
        shell("date > {log}")
        shell("{bcftools} sort  -o {output.vcf} {workdir}/*vcf")
        shell("{bcftools} view -Oz -o {output.vcfgz} {output.vcf}")
