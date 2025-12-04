# Rules for alignment statistics


rule samtools_flagstat:
    """Generate alignment statistics with samtools flagstat."""
    input:
        os.path.join(config["results_dir"], "bam", "{genome}.bam"),
    output:
        os.path.join(config["results_dir"], "stats", "{genome}.flagstat.txt"),
    log:
        os.path.join(config["work_dir"], "logs", "flagstat", "{genome}.log"),
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """


rule samtools_stats:
    """Generate detailed alignment statistics with samtools stats."""
    input:
        os.path.join(config["results_dir"], "bam", "{genome}.bam"),
    output:
        os.path.join(config["results_dir"], "stats", "{genome}.stats.txt"),
    log:
        os.path.join(config["work_dir"], "logs", "stats", "{genome}.log"),
    shell:
        """
        samtools stats {input} > {output} 2> {log}
        """


rule collect_star_logs:
    """Collect STAR alignment logs to results directory."""
    input:
        os.path.join(
            config["work_dir"], "aligned", "{genome}",
            "Log.final.out"
        ),
    output:
        os.path.join(config["results_dir"], "stats", "{genome}.star_log.txt"),
    shell:
        """
        cp {input} {output}
        """


rule multiqc:
    """Aggregate all QC reports with MultiQC."""
    input:
        expand(
            os.path.join(config["results_dir"], "stats", "{genome}.flagstat.txt"),
            genome=get_genomes(),
        ),
        expand(
            os.path.join(config["results_dir"], "stats", "{genome}.star_log.txt"),
            genome=get_genomes(),
        ),
    output:
        os.path.join(config["results_dir"], "multiqc", "multiqc_report.html"),
    log:
        os.path.join(config["work_dir"], "logs", "multiqc.log"),
    params:
        outdir=os.path.join(config["results_dir"], "multiqc"),
        indir=os.path.join(config["results_dir"], "stats"),
    shell:
        """
        multiqc {params.indir} -o {params.outdir} --force 2>&1 | tee {log}
        """
