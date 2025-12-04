# Rules for RNA-seq alignment


rule star_align:
    """Align RNA-seq reads to genome using STAR."""
    input:
        r1=os.path.join(config["work_dir"], "trimmed", "{genome}_R1.fq.gz"),
        r2=os.path.join(config["work_dir"], "trimmed", "{genome}_R2.fq.gz"),
        index=os.path.join(config["work_dir"], "star_index", "{genome}"),
    output:
        bam=os.path.join(
            config["work_dir"], "aligned", "{genome}",
            "Aligned.sortedByCoord.out.bam"
        ),
        log_final=os.path.join(
            config["work_dir"], "aligned", "{genome}",
            "Log.final.out"
        ),
    log:
        os.path.join(config["work_dir"], "logs", "star_align", "{genome}.log"),
    threads: config["star_align"]["threads"]
    resources:
        mem_gb=config["star_align"]["mem_gb"],
    params:
        out_prefix=lambda wildcards, output: output.bam.replace(
            "Aligned.sortedByCoord.out.bam", ""
        ),
        multimap_nmax=config["star_align"]["out_filter_multimap_nmax"],
        mismatch_nmax=config["star_align"]["out_filter_mismatch_nmax"],
    shell:
        """
        STAR --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.out_prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outFilterMultimapNmax {params.multimap_nmax} \
            --outFilterMismatchNmax {params.mismatch_nmax} \
            --outBAMsortingThreadN {threads} \
            2>&1 | tee {log}
        """


rule samtools_index:
    """Index BAM file."""
    input:
        os.path.join(
            config["work_dir"], "aligned", "{genome}",
            "Aligned.sortedByCoord.out.bam"
        ),
    output:
        os.path.join(
            config["work_dir"], "aligned", "{genome}",
            "Aligned.sortedByCoord.out.bam.bai"
        ),
    log:
        os.path.join(config["work_dir"], "logs", "samtools_index", "{genome}.log"),
    shell:
        """
        samtools index {input} 2>&1 | tee {log}
        """


rule copy_final_bam:
    """Copy final BAM and index to results directory."""
    input:
        bam=os.path.join(
            config["work_dir"], "aligned", "{genome}",
            "Aligned.sortedByCoord.out.bam"
        ),
        bai=os.path.join(
            config["work_dir"], "aligned", "{genome}",
            "Aligned.sortedByCoord.out.bam.bai"
        ),
    output:
        bam=os.path.join(config["results_dir"], "bam", "{genome}.bam"),
        bai=os.path.join(config["results_dir"], "bam", "{genome}.bam.bai"),
    shell:
        """
        cp {input.bam} {output.bam}
        cp {input.bai} {output.bai}
        """
