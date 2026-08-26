configfile: "config.yaml"

SAMPLES = ["WT", "AmpB"]

rule all:
    input:
        "results/snpeff/AmpB_unique_nonsyn.tsv",
        expand("results/fastqc/{sample}_{read}_fastqc.html", sample=SAMPLES, read=["1", "2"])

# ---- Session 1: QC, trimming, alignment ----

rule fastqc:
    input:
        r1="data/{sample}_1.fastq.gz",
        r2="data/{sample}_2.fastq.gz"
    output:
        "results/fastqc/{sample}_1_fastqc.html",
        "results/fastqc/{sample}_2_fastqc.html"
    shell:
        "fastqc {input.r1} {input.r2} -o results/fastqc/"

rule trim_galore:
    input:
        r1="data/{sample}_1.fastq.gz",
        r2="data/{sample}_2.fastq.gz"
    output:
        r1="results/trimmed/{sample}_1_val_1.fq.gz",
        r2="results/trimmed/{sample}_2_val_2.fq.gz"
    shell:
        "trim_galore --paired --phred64 -o results/trimmed/ {input.r1} {input.r2}"

rule bowtie2_build:
    input:
        config["reference_fasta"]
    output:
        multiext("results/index/Lmex", ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    shell:
        "bowtie2-build {input} results/index/Lmex"

rule align_sort:
    input:
        r1="results/trimmed/{sample}_1_val_1.fq.gz",
        r2="results/trimmed/{sample}_2_val_2.fq.gz",
        index=rules.bowtie2_build.output
    output:
        bam="results/aligned/{sample}.sorted.bam"
    log:
        "results/aligned/{sample}.bowtie2.log"
    shell:
        # single piped command: align -> BAM -> sort, matching the lab's optional shortcut
        "bowtie2 -x results/index/Lmex --phred64 -1 {input.r1} -2 {input.r2} 2> {log} "
        "| samtools view -b - | samtools sort -o {output.bam} -"

rule index_bam:
    input:
        "results/aligned/{sample}.sorted.bam"
    output:
        "results/aligned/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"

# ---- Session 2: SNP calling ----

rule call_snps:
    input:
        bams=expand("results/aligned/{sample}.sorted.bam", sample=SAMPLES),
        bais=expand("results/aligned/{sample}.sorted.bam.bai", sample=SAMPLES),
        ref=config["reference_fasta"]
    output:
        "results/variants/raw.vcf"
    params:
        ploidy=config["ploidy"]
    shell:
        "bamaddrg -b results/aligned/WT.sorted.bam -s WT "
        "-b results/aligned/AmpB.sorted.bam -s AmpB | "
        "freebayes -f {input.ref} -p {params.ploidy} --stdin > {output}"

# ---- Session 3: filtering, isec, annotation ----

rule filter_snps:
    input:
        "results/variants/raw.vcf"
    output:
        "results/variants/filtered.vcf"
    shell:
        "vcffilter -f 'QUAL > 20' {input} | vcffilter -f 'TYPE = snp' > {output}"

rule split_sample:
    input:
        "results/variants/filtered.vcf"
    output:
        vcf="results/variants/{sample}.vcf.gz",
        idx="results/variants/{sample}.vcf.gz.csi"
    shell:
        "bcftools view -s {wildcards.sample} -c 1 -Oz -o {output.vcf} {input} && "
        "bcftools index {output.vcf}"

rule unique_ampb_snps:
    input:
        wt="results/variants/WT.vcf.gz",
        ampb="results/variants/AmpB.vcf.gz"
    output:
        "results/variants/AmpB_unique.vcf"
    shell:
        # isec writes 4 files into a dir; 0002.vcf = records private to the 2nd file (AmpB)
        "bcftools isec -p results/variants/isec -Ov {input.wt} {input.ampb} && "
        "cp results/variants/isec/0001.vcf {output}"

rule snpeff_annotate:
    input:
        "results/variants/AmpB_unique.vcf"
    output:
        "results/snpeff/AmpB_unique.ann.vcf"
    params:
        config_file=config["snpeff_config"]
    shell:
        "snpEff -Xmx4g -c {params.config_file} -noCheckCds -noCheckProtein -no-intron -no-intergenic "
        "Lmex {input} > {output}"

rule extract_nonsynonymous:
    input:
        "results/snpeff/AmpB_unique.ann.vcf"
    output:
        "results/snpeff/AmpB_unique_nonsyn.tsv"
    shell:
        "cat {input} | vcfEffOnePerLine.pl | "
        "SnpSift extractFields - CHROM POS REF ALT 'ANN[*].IMPACT' 'ANN[*].EFFECT' "
        "'ANN[*].GENE' 'ANN[*].HGVS_C' 'ANN[*].HGVS_P' 'GEN[*].GT' | "
        "grep -E 'HIGH|MODERATE' > {output}"
