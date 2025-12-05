# Rules for alignment statistics


rule samtools_flagstat:
    """Generate alignment statistics with samtools flagstat."""
    input:
        os.path.join(config["results_dir"], "bam", "{sample}.sorted.bam"),
    output:
        os.path.join(config["results_dir"], "stats", "{sample}.flagstat.txt"),
    log:
        os.path.join(config["work_dir"], "logs", "flagstat", "{sample}.log"),
    shell:
        """
        samtools flagstat {input} > {output} 2>&1 | tee {log}
        """


rule count_isoforms:
    """Count isoforms per sample."""
    input:
        gtf=os.path.join(config["results_dir"], "isoforms", "{sample}.isoforms.gtf"),
    output:
        os.path.join(config["results_dir"], "stats", "{sample}.isoform_counts.txt"),
    shell:
        """
        echo "Sample: {wildcards.sample}" > {output}
        echo "Total transcripts: $(grep -c 'transcript' {input.gtf} || echo 0)" >> {output}
        echo "Total exons: $(grep -c 'exon' {input.gtf} || echo 0)" >> {output}
        """


rule multiqc:
    """Generate MultiQC report for all samples."""
    input:
        flagstats=expand(
            os.path.join(config["results_dir"], "stats", "{sample}.flagstat.txt"),
            sample=get_samples(),
        ),
    output:
        os.path.join(config["results_dir"], "multiqc", "multiqc_report.html"),
    log:
        os.path.join(config["work_dir"], "logs", "multiqc.log"),
    params:
        input_dir=os.path.join(config["results_dir"], "stats"),
        output_dir=os.path.join(config["results_dir"], "multiqc"),
    shell:
        """
        multiqc {params.input_dir} \
            --outdir {params.output_dir} \
            --force \
            2>&1 | tee {log}
        """
