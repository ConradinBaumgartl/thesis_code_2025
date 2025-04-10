configfile: "snakemake_scripts/config.yaml"

rule bowtie_align:
    input:
        r1="data/raw_data/{sample}_R1.fastq.gz",
        r2="data/raw_data/{sample}_R2.fastq.gz"
    output:
        cmv_mouse=temp("data/alignment/cmv_mouse/{sample}.sam"),
        cmv_mouse_LOG="data/alignment/cmv_mouse/bowtie2_logs/{sample}.log",
        cmv=temp("data/alignment/cmv/{sample}.sam"),
        cmv_LOG="data/alignment/cmv/bowtie2_logs/{sample}.log"
    params:
        cmv_mouse=config["genomes"]["cmv_mouse"],
        cmv=config["genomes"]["cmv"]
    shell:
        """
        bowtie2 --sensitive-local --no-mixed --no-discordant -p 4 -x {params.cmv_mouse} -1 {input.r1} -2 {input.r2} -S {output.cmv_mouse} &> {output.cmv_mouse_LOG}
        bowtie2 --sensitive-local --no-mixed --no-discordant -p 4 -x {params.cmv} -1 {input.r1} -2 {input.r2} -S {output.cmv} &> {output.cmv_LOG}
        """

rule bam_convert:
    input:
        "data/alignment/{genome}/{sample}.sam"
    output:
        temp("data/alignment/{genome}/{sample}.unsorted.bam")
    shell:
        "samtools view -@ 4 -b {input} > {output}"

rule bam_sort_name:
    input:
        "data/alignment/{genome}/{sample}.unsorted.bam"
    output:
        temp("data/alignment/{genome}/{sample}.n_sorted.bam")
    shell:
        "samtools sort -@ 4 -n {input} > {output}"

rule fixmates:
    input:
        "data/alignment/{genome}/{sample}.n_sorted.bam"
    output:
        temp("data/alignment/{genome}/{sample}.fixmate.tmp")
    shell:
        "samtools fixmate -m {input} -@ 4 {output}"
rule sort_fixmates:
    input:
        "data/alignment/{genome}/{sample}.fixmate.tmp"
    output:
        temp("data/alignment/{genome}/{sample}.fixmate.sorted.tmp")
    shell:
        "samtools sort -@ 4 {input} > {output}"
rule mark_duplicates:
    input:
        "data/alignment/{genome}/{sample}.fixmate.sorted.tmp"
    output:
        "data/alignment/{genome}/{sample}.ALL.bam"
    shell:
        "samtools markdup -@ 4 {input} {output}; samtools index {output}"

rule filter_bam:
    input:
        "data/alignment/{genome}/{sample}.ALL.bam"
    output:
        "data/alignment/{genome}/{sample}.FILT.bam"
    shell:
        "samtools view -b -f 3 -F 3328 -q 30 -@ 4 {input} > {output}; samtools index {output}" # only retain properly paired with alnquality > 30; exclude secondary/supplementary alignment and PCR duplicates

rule size_bam:
    input:
        allreads="data/alignment/{genome}/{sample}.ALL.bam",
        dupfiltered="data/alignment/{genome}/{sample}.FILT.bam"
    output:
        allreads="data/alignment/{genome}/{sample}.SIZE.bam",
        dupfiltered="data/alignment/{genome}/{sample}.FILT.SIZE.bam"
    run:
        shell("alignmentSieve -b {input.allreads} -o {output.allreads} --minFragmentLength 120 --maxFragmentLength 400 -p 4"),
        shell("samtools index {output.allreads}"),
        shell("alignmentSieve -b {input.dupfiltered} -o {output.dupfiltered} --minFragmentLength 120 --maxFragmentLength 400 -p 4"),
        shell("samtools index {output.dupfiltered}")

rule bamqc_ALL:
    input:
        ALL="data/alignment/{genome}/{sample}.ALL.bam",
        FILT="data/alignment/{genome}/{sample}.FILT.bam",
        SIZE="data/alignment/{genome}/{sample}.SIZE.bam",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.bam"
    output:
        ALL_DIR=directory("data/alignment/{genome}/QC/{sample}.ALL/"),
        FILT_DIR=directory("data/alignment/{genome}/QC/{sample}.FILT/"),
        SIZE_DIR=directory("data/alignment/{genome}/QC/{sample}.SIZE/"),
        FILTSIZE_DIR=directory("data/alignment/{genome}/QC/{sample}.FILT.SIZE/"),
        ALL="data/alignment/{genome}/QC/{sample}.ALL/genome_results.txt",
        FILT="data/alignment/{genome}/QC/{sample}.FILT/genome_results.txt",
        SIZE="data/alignment/{genome}/QC/{sample}.SIZE/genome_results.txt",
        FILTSIZE="data/alignment/{genome}/QC/{sample}.FILT.SIZE/genome_results.txt"
    run:
        shell("qualimap bamqc -bam {input.ALL} -outdir {output.ALL_DIR}"),
        shell("qualimap bamqc -bam {input.FILT} -outdir {output.FILT_DIR}"),
        shell("qualimap bamqc -bam {input.SIZE} -outdir {output.SIZE_DIR}"),
        shell("qualimap bamqc -bam {input.FILTSIZE} -outdir {output.FILTSIZE_DIR}")

rule bam_coverage_ALL:
    input:
        ALL="data/alignment/{genome}/{sample}.ALL.bam",
        FILT="data/alignment/{genome}/{sample}.FILT.bam",
        SIZE="data/alignment/{genome}/{sample}.SIZE.bam",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.bam"
    output:
        ALL="data/alignment/{genome}/{sample}.ALL.bw",
        FILT="data/alignment/{genome}/{sample}.FILT.bw",
        SIZE="data/alignment/{genome}/{sample}.SIZE.bw",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.bw"
    run:
        shell("bamCoverage --normalizeUsing CPM --extendReads -bs 1 -p 4 -b {input.ALL} -o {output.ALL}"),
        shell("bamCoverage --normalizeUsing CPM --extendReads -bs 1 -p 4 -b {input.FILT} -o {output.FILT}"),
        shell("bamCoverage --normalizeUsing CPM --extendReads -bs 1 -p 4 -b {input.SIZE} -o {output.SIZE}"),
        shell("bamCoverage --normalizeUsing CPM --extendReads -bs 1 -p 4 -b {input.FILTSIZE} -o {output.FILTSIZE}")

rule bam_coverage_rawnumbers:
    input:
        ALL="data/alignment/{genome}/{sample}.ALL.bam",
        FILT="data/alignment/{genome}/{sample}.FILT.bam",
        SIZE="data/alignment/{genome}/{sample}.SIZE.bam",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.bam"
    output:
        ALL="data/alignment/{genome}/{sample}.ALL.nreads",
        FILT="data/alignment/{genome}/{sample}.FILT.nreads",
        SIZE="data/alignment/{genome}/{sample}.SIZE.nreads",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.nreads"
    run:
        shell("samtools idxstats {input.ALL} > {output.ALL}"),
        shell("samtools idxstats {input.SIZE} > {output.SIZE}"),
        shell("samtools idxstats {input.FILT} > {output.FILT}"),
        shell("samtools idxstats {input.FILTSIZE} > {output.FILTSIZE}")


rule insert_size_metrics:
    input:
        ALL="data/alignment/{genome}/{sample}.ALL.bam",
        FILT="data/alignment/{genome}/{sample}.FILT.bam",
        SIZE="data/alignment/{genome}/{sample}.SIZE.bam",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.bam"
    output:
        ALLH="data/alignment/{genome}/{sample}.ALL.H",
        FILTH="data/alignment/{genome}/{sample}.FILT.H",
        SIZEH="data/alignment/{genome}/{sample}.SIZE.H",
        FILTSIZEH="data/alignment/{genome}/{sample}.FILT.SIZE.H",
        ALLhist="data/alignment/{genome}/{sample}.ALL.histogram",
        FILThist="data/alignment/{genome}/{sample}.FILT.histogram",
        SIZEhist="data/alignment/{genome}/{sample}.SIZE.histogram",
        FILTSIZEhist="data/alignment/{genome}/{sample}.FILT.SIZE.histogram"
    run:
        shell("picard CollectInsertSizeMetrics --INCLUDE_DUPLICATES true -I {input.ALL} -O {output.ALLhist} -H {output.ALLH}"),
        shell("picard CollectInsertSizeMetrics --INCLUDE_DUPLICATES true -I {input.FILT} -O {output.FILThist} -H {output.FILTH}"),
        shell("picard CollectInsertSizeMetrics --INCLUDE_DUPLICATES true -I {input.SIZE} -O {output.SIZEhist} -H {output.SIZEH}"),
        shell("picard CollectInsertSizeMetrics --INCLUDE_DUPLICATES true -I {input.FILTSIZE} -O {output.FILTSIZEhist} -H {output.FILTSIZEH}")



#rule gather_info:
#    input:
#        ALL="data/alignment/{genome}/QC/{sample}/genome_results.txt",
#        FILT="data/alignment/{genome}/QC/{sample}.FILT/genome_results.txt",
#        SIZE="data/alignment/{genome}/QC/{sample}.SIZE/genome_results.txt",
#        FILTSIZE="data/alignment/{genome}/QC/{sample}.FILT.SIZE/genome_results.txt"
#    output:
#
#    run:
#        shell("grep ’number of reads’ {input} | awk 'BEGIN {FS = '='; OFS = "\\t"} {gsub("\\s", "_", $1); gsub(",", "", $2); print substr($1, 6, length($1)-6) , $2}' >> data/alignment/QC/${sample}.info")
#
#  grep "number of mapped reads" data/alignment/QC/${sample}/genome_results.txt | awk 'BEGIN {FS = "="; OFS = "\t"} {gsub("\\s", "_", $1); gsub(",", "", $2); print substr($1, 6, length($1)-6) , substr($2, 1, #length($2)-8)}' >> data/alignment/QC/${sample}.info
#  grep "number of duplicated" data/alignment/QC/${sample}/genome_results.txt | awk 'BEGIN {FS = "="; OFS = "\t"} {gsub("\\s", "_", $1); gsub(",", "", $2); print substr($1, 6, length($1)-16), $2}' >> data/#alignment/QC/${sample}.info
#  grep "mean mapping quality" data/alignment/QC/${sample}/genome_results.txt | awk 'BEGIN {FS = "="; OFS = "\t"} {gsub("\\s", "_", $1); gsub(",", "", $2); print substr($1, 6, length($1)-6), $2}' >> data/#alignment/QC/${sample}.info
#  grep "pPl04" data/alignment/${sample}.nreads | awk 'BEGIN {FS = "\t"; OFS = "\t"} {print "number_of_AAV_mapped_reads", $3}' >> data/alignment/QC/${sample}.info
