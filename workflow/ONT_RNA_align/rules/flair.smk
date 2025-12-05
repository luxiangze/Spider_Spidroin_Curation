# Rules for FLAIR workflow
# Using flair transcriptome (combines correct + collapse) with short-read junction support


rule gff_to_gtf:
    """
    Convert GFF3 annotation to GTF format for FLAIR compatibility.
    FLAIR requires GTF format with transcript_id attribute.
    """
    input:
        gff=get_genome_gff,
    output:
        gtf=os.path.join(config["work_dir"], "annotation", "{sample}.gtf"),
    log:
        os.path.join(config["work_dir"], "logs", "gff_to_gtf", "{sample}.log"),
    shell:
        """
        gffread {input.gff} -T -o {output.gtf} 2> {log}
        """


rule minimap2_align:
    """
    Step 1: Align ONT reads to genome using minimap2.
    Recommended alignment for FLAIR transcriptome input.
    """
    input:
        reads=get_fastq,
        genome=get_genome_fasta,
    output:
        bam=os.path.join(config["work_dir"], "aligned", "{sample}.bam"),
    log:
        os.path.join(config["work_dir"], "logs", "minimap2", "{sample}.log"),
    threads: config["minimap2"]["threads"]
    resources:
        mem_gb=config["minimap2"]["mem_gb"],
    params:
        extra=config["minimap2"].get("extra", ""),
    shell:
        """
        minimap2 -ax splice {params.extra} \
            -t {threads} \
            --MD --secondary=no \
            {input.genome} {input.reads} 2> {log} | \
            samtools view -hb - | \
            samtools sort -@ 4 -o {output.bam} -
        samtools index -c {output.bam}
        """


rule flair_transcriptome:
    """
    Step 2: Build transcriptome with FLAIR (combines correct + collapse).
    Uses short-read junction support from BGI STAR alignment.
    """
    input:
        bam=os.path.join(config["work_dir"], "aligned", "{sample}.bam"),
        genome=get_genome_fasta,
        gtf=os.path.join(config["work_dir"], "annotation", "{sample}.gtf"),
        junction_tab=get_star_sj_tab if config.get("use_short_read_junctions", False) else [],
    output:
        gtf=os.path.join(config["work_dir"], "transcriptome", "{sample}.isoforms.gtf"),
        fa=os.path.join(config["work_dir"], "transcriptome", "{sample}.isoforms.fa"),
        bed=os.path.join(config["work_dir"], "transcriptome", "{sample}.isoforms.bed"),
    log:
        os.path.join(config["work_dir"], "logs", "flair_transcriptome", "{sample}.log"),
    threads: config["flair_transcriptome"]["threads"]
    resources:
        mem_gb=config["flair_transcriptome"]["mem_gb"],
    params:
        out_prefix=lambda wildcards: os.path.join(
            config["work_dir"], "transcriptome", wildcards.sample
        ),
        stringent=lambda wildcards: "--stringent" if config["flair_transcriptome"]["stringent"] else "",
        min_support=config["flair_transcriptome"]["min_support"],
        junction_opt=lambda wildcards, input: f"--junction_tab {input.junction_tab}" if input.junction_tab else "",
    shell:
        """
        flair transcriptome \
            -g {input.genome} \
            -f {input.gtf} \
            -b {input.bam} \
            -o {params.out_prefix} \
            -t {threads} \
            --support {params.min_support} \
            {params.junction_opt} \
            {params.stringent} \
            2>&1 | tee {log}
        """


rule copy_isoforms:
    """Copy isoform results to results directory."""
    input:
        gtf=os.path.join(config["work_dir"], "transcriptome", "{sample}.isoforms.gtf"),
        fa=os.path.join(config["work_dir"], "transcriptome", "{sample}.isoforms.fa"),
        bed=os.path.join(config["work_dir"], "transcriptome", "{sample}.isoforms.bed"),
    output:
        gtf=os.path.join(config["results_dir"], "isoforms", "{sample}.isoforms.gtf"),
        fa=os.path.join(config["results_dir"], "isoforms", "{sample}.isoforms.fa"),
        bed=os.path.join(config["results_dir"], "isoforms", "{sample}.isoforms.bed"),
    shell:
        """
        cp {input.gtf} {output.gtf}
        cp {input.fa} {output.fa}
        cp {input.bed} {output.bed}
        """


rule copy_bam:
    """Copy aligned BAM to results directory."""
    input:
        os.path.join(config["work_dir"], "aligned", "{sample}.bam"),
    output:
        os.path.join(config["results_dir"], "bam", "{sample}.bam"),
    shell:
        """
        cp {input} {output}
        """
