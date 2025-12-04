# Rules for quality control


rule fastqc_raw:
    """Run FastQC on raw reads."""
    input:
        r1=get_fastq_r1,
        r2=get_fastq_r2,
    output:
        html1=os.path.join(config["work_dir"], "fastqc_raw", "{genome}_R1_fastqc.html"),
        html2=os.path.join(config["work_dir"], "fastqc_raw", "{genome}_R2_fastqc.html"),
        zip1=os.path.join(config["work_dir"], "fastqc_raw", "{genome}_R1_fastqc.zip"),
        zip2=os.path.join(config["work_dir"], "fastqc_raw", "{genome}_R2_fastqc.zip"),
    log:
        os.path.join(config["work_dir"], "logs", "fastqc_raw", "{genome}.log"),
    threads: 2
    params:
        outdir=os.path.join(config["work_dir"], "fastqc_raw"),
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} 2>&1 | tee {log}
        # Rename to consistent naming
        mv {params.outdir}/$(basename {input.r1} .fq.gz)_fastqc.html {output.html1} || true
        mv {params.outdir}/$(basename {input.r2} .fq.gz)_fastqc.html {output.html2} || true
        mv {params.outdir}/$(basename {input.r1} .fq.gz)_fastqc.zip {output.zip1} || true
        mv {params.outdir}/$(basename {input.r2} .fq.gz)_fastqc.zip {output.zip2} || true
        """


rule fastp:
    """Quality control and trimming with fastp."""
    input:
        r1=get_fastq_r1,
        r2=get_fastq_r2,
    output:
        r1=os.path.join(config["work_dir"], "trimmed", "{genome}_R1.fq.gz"),
        r2=os.path.join(config["work_dir"], "trimmed", "{genome}_R2.fq.gz"),
        html=os.path.join(config["work_dir"], "fastp", "{genome}.fastp.html"),
        json=os.path.join(config["work_dir"], "fastp", "{genome}.fastp.json"),
    log:
        os.path.join(config["work_dir"], "logs", "fastp", "{genome}.log"),
    threads: config["fastp"]["threads"]
    params:
        extra=config["fastp"].get("extra", ""),
    shell:
        """
        fastp -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -h {output.html} -j {output.json} \
            -w {threads} \
            {params.extra} \
            2>&1 | tee {log}
        """


rule fastqc_trimmed:
    """Run FastQC on trimmed reads."""
    input:
        r1=os.path.join(config["work_dir"], "trimmed", "{genome}_R1.fq.gz"),
        r2=os.path.join(config["work_dir"], "trimmed", "{genome}_R2.fq.gz"),
    output:
        html1=os.path.join(config["work_dir"], "fastqc_trimmed", "{genome}_R1_fastqc.html"),
        html2=os.path.join(config["work_dir"], "fastqc_trimmed", "{genome}_R2_fastqc.html"),
        zip1=os.path.join(config["work_dir"], "fastqc_trimmed", "{genome}_R1_fastqc.zip"),
        zip2=os.path.join(config["work_dir"], "fastqc_trimmed", "{genome}_R2_fastqc.zip"),
    log:
        os.path.join(config["work_dir"], "logs", "fastqc_trimmed", "{genome}.log"),
    threads: 2
    params:
        outdir=os.path.join(config["work_dir"], "fastqc_trimmed"),
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} 2>&1 | tee {log}
        mv {params.outdir}/$(basename {input.r1} .fq.gz)_fastqc.html {output.html1} || true
        mv {params.outdir}/$(basename {input.r2} .fq.gz)_fastqc.html {output.html2} || true
        mv {params.outdir}/$(basename {input.r1} .fq.gz)_fastqc.zip {output.zip1} || true
        mv {params.outdir}/$(basename {input.r2} .fq.gz)_fastqc.zip {output.zip2} || true
        """


rule multiqc_raw:
    """MultiQC report for raw reads QC."""
    input:
        expand(
            os.path.join(config["work_dir"], "fastqc_raw", "{genome}_{read}_fastqc.zip"),
            genome=get_genomes(),
            read=["R1", "R2"],
        ),
    output:
        os.path.join(config["results_dir"], "multiqc", "raw_reads_multiqc.html"),
    log:
        os.path.join(config["work_dir"], "logs", "multiqc_raw.log"),
    params:
        indir=os.path.join(config["work_dir"], "fastqc_raw"),
        outdir=os.path.join(config["results_dir"], "multiqc"),
    shell:
        """
        multiqc {params.indir} \
            -o {params.outdir} \
            -n raw_reads_multiqc.html \
            --force 2>&1 | tee {log}
        """


rule multiqc_trimmed:
    """MultiQC report for trimmed reads QC."""
    input:
        fastqc=expand(
            os.path.join(config["work_dir"], "fastqc_trimmed", "{genome}_{read}_fastqc.zip"),
            genome=get_genomes(),
            read=["R1", "R2"],
        ),
        fastp=expand(
            os.path.join(config["work_dir"], "fastp", "{genome}.fastp.json"),
            genome=get_genomes(),
        ),
    output:
        os.path.join(config["results_dir"], "multiqc", "trimmed_reads_multiqc.html"),
    log:
        os.path.join(config["work_dir"], "logs", "multiqc_trimmed.log"),
    params:
        fastqc_dir=os.path.join(config["work_dir"], "fastqc_trimmed"),
        fastp_dir=os.path.join(config["work_dir"], "fastp"),
        outdir=os.path.join(config["results_dir"], "multiqc"),
    shell:
        """
        multiqc {params.fastqc_dir} {params.fastp_dir} \
            -o {params.outdir} \
            -n trimmed_reads_multiqc.html \
            --force 2>&1 | tee {log}
        """
