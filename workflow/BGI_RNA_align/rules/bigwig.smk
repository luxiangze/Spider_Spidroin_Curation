# Rules for generating bigWig files


rule bam_to_bigwig:
    """Convert BAM to bigWig for genome browser visualization."""
    input:
        bam=os.path.join(config["results_dir"], "bam", "{genome}.bam"),
        csi=os.path.join(config["results_dir"], "bam", "{genome}.bam.csi"),
    output:
        os.path.join(config["results_dir"], "bigwig", "{genome}.bw"),
    log:
        os.path.join(config["work_dir"], "logs", "bigwig", "{genome}.log"),
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


rule bam_to_bigwig_stranded:
    """Generate strand-specific bigWig files."""
    input:
        bam=os.path.join(config["results_dir"], "bam", "{genome}.bam"),
        csi=os.path.join(config["results_dir"], "bam", "{genome}.bam.csi"),
    output:
        fwd=os.path.join(config["results_dir"], "bigwig", "{genome}.fwd.bw"),
        rev=os.path.join(config["results_dir"], "bigwig", "{genome}.rev.bw"),
    log:
        os.path.join(config["work_dir"], "logs", "bigwig_stranded", "{genome}.log"),
    threads: config["bigwig"]["threads"]
    params:
        normalize=config["bigwig"]["normalize_using"],
        bin_size=config["bigwig"]["bin_size"],
    shell:
        """
        # Forward strand
        bamCoverage -b {input.bam} -o {output.fwd} \
            --binSize {params.bin_size} \
            --normalizeUsing {params.normalize} \
            --filterRNAstrand forward \
            --numberOfProcessors {threads} \
            2>&1 | tee {log}

        # Reverse strand
        bamCoverage -b {input.bam} -o {output.rev} \
            --binSize {params.bin_size} \
            --normalizeUsing {params.normalize} \
            --filterRNAstrand reverse \
            --numberOfProcessors {threads} \
            2>&1 >> {log}
        """
