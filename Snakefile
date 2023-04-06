ALL_STRAINS = ['DL238', 'DRR142768', 'ECA36', 'ECA396', 'EG4725', 'JU1400', 'JU2526', 'JU2600', 'JU310', 'MY2147', 'MY2693', 'NIC2', 'NIC526', 'QX1794', 'XZ1516']
REFERENCE = "0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa" # C. elegans reference genome
NGMLRDICT40X = {"DL238": "0.24", "DRR142768": "0.52", "ECA36": "0.19", "ECA396": "0.30", "EG4725": "0.18", "JU1400": "0.34", "JU2526": "0.33", "JU2600": "0.22", "JU310": "0.25", "MY2147":"0.21", "MY2693": "0.31", "NIC2": "0.33", "NIC526": "0.27", "QX1794": "0.32","XZ1516": "0.54"}

rule all:
	input:
		expand("1_alignments/ngmlr/{strain}/{strain}_picard_sorted.{extension}", strain=ALL_STRAINS, extension=["bam","bai"]),

# Align FASTQ files using NGMLR
#rule ngmlr:
#	input:
#		ancient("0_input/fastq/{strain}_all_reads.fastq")
#	output:
#		temp("1_alignments/ngmlr/{strain}/{strain}_aln.sam")
#	conda:  "yaml/ngmlr.yaml"
#	threads: 8
#	resources:
#		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
#		time="60:00:00"
#	shell:
#		"ngmlr -t {threads} -r {REFERENCE} -q {input} -o {output}"

rule picard_sort:
	input:
			"1_alignments/ngmlr/{strain}/{strain}_sorted.bam"
	        #"1_alignments/ngmlr/{strain}/{strain}_aln.sam"
	output:
	        bamfile="1_alignments/ngmlr/{strain}/{strain}_picard_sorted.bam",
	        bamindex="1_alignments/ngmlr/{strain}/{strain}_picard_sorted.bai"
	params:
	        picard_cmd=r"""java "-Xmx60g" -jar /home/kyle.lesack1/miniconda3/envs/picardtools/share/picard-2.27.5-0/picard.jar SortSam """,
			max_records="25000"
	conda:  "yaml/picardtools.v2.27.5.yaml"
	threads: 8
	resources:
	        mem_mb=lambda _, attempt: 60000 + ((attempt - 1) * 10000),
	        time_hms="02:00:00"
	shell:
	        "{params.picard_cmd} -I {input} -O {output.bamfile} -SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT --CREATE_INDEX --MAX_RECORDS_IN_RAM {params.max_records}"

# Rename bam index file to bam.bai in order to work with SVIM
rule rename_bam_index:
	input:
	        "1_alignments/ngmlr/{strain}/{strain}_picard_sorted.bai"
	output:
	        "1_alignments/ngmlr/{strain}/{strain}_picard_sorted.bam.bai"
	conda:  "yaml/rename.yaml"
	threads: 1
	resources:
		mem_mb=100,
		time_hms="00:05:00"
	shell:
		"rename 's/\.bai/.bam.bai/' {input}"

# Call SVs with SVIM
rule svim:
	input:
	        bamfile="1_alignments/ngmlr/{strain}/{strain}_picard_sorted.bam",
	        bamindex="1_alignments/ngmlr/{strain}/{strain}_picard_sorted.bai"
	output:
		"2_variant_calls/ngmlr/svim/{strain}/variants.vcf"
	params:
		outdir="2_variant_calls/ngmlr/svim/{strain}/variants.vcf"
	conda:  "yaml/svim2.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 20000 + ((attempt - 1) * 10000),
		time_hms="02:00:00"
	shell:
		"svim alignment {params.outdir} {input.bamfile} {REFERENCE}"
