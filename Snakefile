ALL_STRAINS = ["JU1400", "NIC2", "JU2526", "XZ1516", "MY2693", "QX1794", "NIC526", "DL238", "ECA396","JU2600","ECA36","EG4725","MY2147","JU310"]
REFERENCE = "0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa" # C. elegans reference genome
SAM_ALIGNERS = ["ngmlr","minimap2"] # Aligners that output sam files

rule all:
	input:
		expand("1_alignments/minimap2/{strain}/{strain}_sorted.{extension}", strain=ALL_STRAINS, extension=["bam","bam.bai"]),

# Align FASTQ files using Minimap2
rule minimap2_old:
	input:
		ancient("0_input/fastq/{strain}_all_reads.fastq")
	output:
		temp("2_alignments/full_depth/minimap2/{strain}/{readorder}/{strain}_{readorder}.sam")
	conda:  "yaml/pbhoney.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="4-0:00:00"
	shell:
		"minimap2 -t {threads} -ax map-pb {REFERENCE} {input} > {output}"

rule minimap2_index:
    input:
        target="0_data/reference/c_elegans.PRJNA13758.WS263.genomic.fa"
    output:
        "0_data/reference/c_elegans.PRJNA13758.WS263.genomic.fa.mmi"
    log:
        "logs/minimap2_index/{strain}.log"
    params:
        extra=""  # optional additional args
    threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="02:00:00"
    wrapper:
        "v1.25.0/bio/minimap2/index"

rule minimap2_sam:
    input:
        target="0_data/reference/c_elegans.PRJNA13758.WS263.genomic.fa.mmi",  # can be either genome index or genome fasta
        query="0_input/fastq/{strain}_all_reads.fastq"
    output:
        temp("1_alignments/minimap2/{strain}/{strain}_aln.sam")
    log:
        "logs/minimap2/{strain}.log",
    params:
        extra="-ax map-pb",  # optional
        sorting="none",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra="",  # optional: extra arguments for samtools/picard
    threads: 8
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time="4-0:00:00"
    wrapper:
        "v1.25.0/bio/minimap2/aligner"

rule picard_sort:
        input:
                "1_alignments/minimap2/{strain}/{strain}_aln.sam"
        output:
                "1_alignments/minimap2/{strain}/{strain}_sorted.bam",
                "1_alignments/minimap2/{strain}/{strain}_sorted.bam.bai"
        params:
                picard_cmd=r"""java "-Xmx70g" -jar /home/kyle.lesack1/miniconda3/envs/picardtools/share/picard-2.27.5-0/picard.jar SortSam """,
                outbam = "2_alignments.new/full_depth/{aligner}/{strain}/shuffled/{strain}_shuffled_picardsorted.bam"
        conda:  "yaml/picardtools.v2.27.5.yaml"
        threads: 8
        resources:
                mem_mb=lambda _, attempt: 40000 + ((attempt - 1) * 10000),
                time_hms="04:00:00"
        shell:
                "{params.picard_cmd} -I {input} -O {params.outbam} -SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT --CREATE_INDEX --MAX_RECORDS_IN_RAM 250000"
