workdir: "/scif/data"

rule bwa_index:
    input:
        "index/{ref}.fa"
    output:
        "index/{ref}.fa.bwt"
    shell:
        "scif run bwa index {input}"

rule bwa_map:
    input:
        "index/ref.fa",
        "raw_reads/{sample}.fq"
    output:
        "mapped_reads/{sample}.sam"
    shell:
        "scif run bwa mem {input} > {output}"


rule samtools_sam_to_bam:
    input:
        "mapped_reads/{sample}.sam"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "scif run samtools view -Sb {input} > {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "scif run samtools sort -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "scif run samtools index {input}"
