rule bam_coverage_bedgraph:
    input:
        ALL="data/alignment/{genome}/{sample}.ALL.bam",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.bam"
    output:
        ALL="data/alignment/{genome}/{sample}.ALL.bedgraph",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.bedgraph"
    run:
        shell("bamCoverage -of bedgraph --extendReads -bs 1 -p 4 --samFlagInclude 67 -b {input.ALL} -o {output.ALL}"),
        shell("bamCoverage -of bedgraph --extendReads -bs 1 -p 4 --samFlagInclude 67 -b {input.FILTSIZE} -o {output.FILTSIZE}")

rule SEACR:
    input:
        ALL="data/alignment/{genome}/{sample}.ALL.bedgraph",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.bedgraph"
    params:
        ALL="data/alignment/{genome}/{sample}.ALL.SEACR",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.SEACR"
    output:
        ALL="data/alignment/{genome}/{sample}.ALL.SEACR.stringent.bed",
        FILTSIZE="data/alignment/{genome}/{sample}.FILT.SIZE.SEACR.stringent.bed"
    run:
        shell("bash ~/applications/SEACR/SEACR_1.3.sh {input.ALL} 0.01 norm stringent {params.ALL}"),
        shell("bash ~/applications/SEACR/SEACR_1.3.sh {input.FILTSIZE} 0.01 norm stringent {params.FILTSIZE}")
