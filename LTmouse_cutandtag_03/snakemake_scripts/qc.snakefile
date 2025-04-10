rule fastqc:
    input:
        r1=expand("data/raw_data/{samples}_R1.fastq.gz", samples = samples),
        r2=expand("data/raw_data/{samples}_R2.fastq.gz", samples = samples)
    output:
        expand("data/raw_data/QC/{samples}_R1_fastqc.html", samples=samples),
        expand("data/raw_data/QC/{samples}_R2_fastqc.html", samples=samples)
    params:
        outdir="data/raw_data/QC",
        threads=threads
    benchmark:
        "data/raw_data/QC/logs/fastqc_benchmark.txt"
    shell:
        "fastqc {input} -o {params.outdir} -t {params.threads}"


rule multiqc:
    input:
        expand("data/raw_data/QC/{samples}_R1_fastqc.html", samples=samples),
        expand("data/raw_data/QC/{samples}_R2_fastqc.html", samples=samples),
        "data/raw_data/stats"
    output:
        "data/raw_data/QC/multiqc/multiqc_report.html"
    params:
        indir="data/raw_data/QC/.",
        outdir="data/raw_data/QC/multiqc/"
    benchmark:
        "data/raw_data/QC/logs/multiqc_benchmark.log"
    shell:
        "multiqc {params.indir} -s -o {params.outdir}"

rule seqkit_stats:
    input:
        r1=expand("data/raw_data/{samples}_R1.fastq.gz", samples = samples),
        r2=expand("data/raw_data/{samples}_R2.fastq.gz", samples = samples)
    output:
        "data/raw_data/stats"
    shell:
        "seqkit stats data/raw_data/*.fastq.gz > data/raw_data/stats"
