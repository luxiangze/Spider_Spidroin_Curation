# Rules for building genome indices


rule star_index:
    """Build STAR index for a genome."""
    input:
        fasta=get_genome_fasta,
        gff=get_genome_gff,
    output:
        directory(os.path.join(config["work_dir"], "star_index", "{genome}")),
    log:
        os.path.join(config["work_dir"], "logs", "star_index", "{genome}.log"),
    threads: config["star_index"]["threads"]
    resources:
        mem_gb=config["star_index"]["mem_gb"],
    params:
        sjdb_overhang=config["star_index"]["sjdb_overhang"],
        sa_index_n_bases=config["star_index"]["genome_sa_index_n_bases"],
    shell:
        """
        mkdir -p {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gff} \
            --sjdbGTFtagExonParentTranscript Parent \
            --sjdbOverhang {params.sjdb_overhang} \
            --genomeSAindexNbases {params.sa_index_n_bases} \
            2>&1 | tee {log}
        """


rule samtools_faidx:
    """Create fasta index with samtools."""
    input:
        get_genome_fasta,
    output:
        os.path.join(config["work_dir"], "faidx", "{genome}.fa.fai"),
    log:
        os.path.join(config["work_dir"], "logs", "samtools_faidx", "{genome}.log"),
    shell:
        """
        samtools faidx {input} -o {output} 2>&1 | tee {log}
        """
