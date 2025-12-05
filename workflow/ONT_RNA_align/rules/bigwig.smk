# Rules for generating bigWig files from ONT alignments


rule samtools_sort:
    """Sort BAM file for bigWig generation."""
    input:
        os.path.join(config["results_dir"], "bam", "{sample}.bam"),
    output:
        os.path.join(config["results_dir"], "bam", "{sample}.sorted.bam"),
    log:
        os.path.join(config["work_dir"], "logs", "samtools_sort", "{sample}.log"),
    threads: 8
    shell:
        """
        samtools sort -@ {threads} -o {output} {input} 2>&1 | tee {log}
        """


rule samtools_index:
    """Index sorted BAM file using CSI index for large genomes."""
    input:
        os.path.join(config["results_dir"], "bam", "{sample}.sorted.bam"),
    output:
        os.path.join(config["results_dir"], "bam", "{sample}.sorted.bam.csi"),
    log:
        os.path.join(config["work_dir"], "logs", "samtools_index", "{sample}.log"),
    shell:
        """
        samtools index -c {input} 2>&1 | tee {log}
        """


rule bam_to_bigwig:
    """Convert BAM to bigWig for genome browser visualization."""
    input:
        bam=os.path.join(config["results_dir"], "bam", "{sample}.sorted.bam"),
        csi=os.path.join(config["results_dir"], "bam", "{sample}.sorted.bam.csi"),
    output:
        os.path.join(config["results_dir"], "bigwig", "{sample}.bw"),
    log:
        os.path.join(config["work_dir"], "logs", "bigwig", "{sample}.log"),
    threads: config["bigwig"]["threads"]
    params:
        normalize=config["bigwig"]["normalize_using"],
        bin_size=config["bigwig"]["bin_size"],
    shell:
        """
        bamCoverage -b {input.bam} -o {output} \
            --binSize {params.bin_size} \
            --normalizeUsing {params.normalize} \
            --numberOfProcessors {threads} \
            2>&1 | tee {log}
        """
