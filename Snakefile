ALL_STRAINS = ['DL238', 'DRR142768', 'ECA36', 'ECA396', 'EG4725', 'JU1400', 'JU2526', 'JU2600', 'JU310', 'MY2147', 'MY2693', 'NIC2', 'NIC526', 'QX1794', 'XZ1516']
REFERENCE = "0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa" # C. elegans reference genome
NGMLRDICT40X = {"DL238": "0.24", "DRR142768": "0.52", "ECA36": "0.19", "ECA396": "0.30", "EG4725": "0.18", "JU1400": "0.34", "JU2526": "0.33", "JU2600": "0.22", "JU310": "0.25", "MY2147":"0.21", "MY2693": "0.31", "NIC2": "0.33", "NIC526": "0.27", "QX1794": "0.32","XZ1516": "0.54"}
ALIGNMENT_DIR = ["ngmlr", "ngmlr.subsampled_40X"]

rule all:
	input:
		expand("1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.{extension}", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS, extension=["bam","bam.bai"]),
		#expand("1_alignments/ngmlr.subsampled_40X/{strain}/{strain}_picard_sorted.bam", strain=ALL_STRAINS),
		#expand("1_alignments/ngmlr.subsampled.picard_sorted/{strain}/{strain}_picard_sorted.bam", strain=ALL_STRAINS), # use diff on these and subsampled. Remove this if same.
		expand("2_variant_calls/{alignment_dir}/svim/{strain}/variants.vcf", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS),
		#expand("2_variant_calls/ngmlr.subsampled_40X/svim/{strain}/variants.vcf", strain=ALL_STRAINS),
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
			"1_alignments/{alignment_dir}/{strain}/{strain}_sorted.bam"
	        #"1_alignments/ngmlr/{strain}/{strain}_aln.sam"
	output:
	        bamfile="1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.bam",
	        bamindex="1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.bai"
	params:
	        picard_cmd=r"""java "-Xmx60g" -jar /home/kyle.lesack1/miniconda3/envs/picardtools/share/picard-2.27.5-0/picard.jar SortSam """,
			max_records="25000",
	conda:  "yaml/picardtools.v2.27.5.yaml"
	#conda:  "yaml/picardtools.v2.27.5_rename" #rename 's/\.bai/.bam.bai/' {params.tempindex}
	threads: 8
	resources:
	        mem_mb=lambda _, attempt: 60000 + ((attempt - 1) * 10000),
	        time_hms="02:00:00"
	shell:
	        "{params.picard_cmd} -I {input} -O {output.bamfile} -SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT --CREATE_INDEX --MAX_RECORDS_IN_RAM {params.max_records}"

# I can't seem to get the rename command to work with the picard sort rule. Might be a Slurm issue
rule copy_bam_index:
	input:
			"1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.bai"

	output:
	        "1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.bam.bai"
	threads: 1
	resources:
	        mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
	        time_hms="00:05:00"
	shell:
	        "cp {input} {output}"

# Call SVs with SVIM
rule svim:
	input:
	        bamfile="1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.bam",
	        bamindex="1_alignments/ngmlr/{strain}/{strain}_picard_sorted.bam.bai"
	output:
		"2_variant_calls/{alignment_dir}/svim/{strain}/variants.vcf"
	params:
		outdir="2_variant_calls/{alignment_dir}/svim/{strain}/",
		MIN_SV_SIZE="100",
		MINIMUM_SCORE="15"
	conda:  "yaml/svim2.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 20000 + ((attempt - 1) * 10000),
		time_hms="02:00:00"
	shell:
		"svim alignment --min_sv_size {params.MIN_SV_SIZE} --minimum_score {params.MINIMUM_SCORE} {params.outdir} {input.bamfile} {REFERENCE}"

# Call SVs with SVIM
#rule svim_subsampled:
#	input:
#	        bamfile="1_alignments/ngmlr.subsampled_40X/{strain}/{strain}_picard_sorted.bam",
#	        bamindex="1_alignments/ngmlr.subsampled_40X/{strain}/{strain}_picard_sorted.bam.bai"
#	output:
#		"2_variant_calls/ngmlr.subsampled_40X/svim/{strain}/variants.vcf"
#	params:
#		outdir="2_variant_calls/ngmlr.subsampled_40X/svim/{strain}/",
#		MIN_SV_SIZE="100",
#		MINIMUM_SCORE="15"
#	conda:  "yaml/svim2.yaml"
#	threads: 8
#	resources:
#		mem_mb=lambda _, attempt: 20000 + ((attempt - 1) * 10000),
#		time_hms="02:00:00"
#	shell:
#		"svim alignment --min_sv_size {params.MIN_SV_SIZE} --minimum_score {params.MINIMUM_SCORE} {params.outdir} {input.bamfile} {REFERENCE}"

rule subsample_ngmlr_40x:
	input:
		"1_alignments/ngmlr/{strain}/{strain}_picard_sorted.bam"
	output:
		"1_alignments/ngmlr.subsampled_40X/{strain}/{strain}_sorted.bam"
	params:
		value=lambda wcs: NGMLRDICT40X[wcs.strain]
	conda:  "yaml/samtools_1.9.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 10000 + ((attempt - 1) * 10000),
		time_hms="08:00:00"
	shell:
		"""samtools view -@ {threads} -b -s {params.value} {input} > {output}"""

# Check if the subsample_ngmlr_40x BAM files remain sorted
#rule picard_sort_subsampled:
#	input:
#			"1_alignments/ngmlr.subsampled_40X/{strain}/{strain}_picard_sorted.bam"
#	output:
#	        bamfile="1_alignments/ngmlr.subsampled.picard_sorted/{strain}/{strain}_picard_sorted.bam",
#	params:
#	        picard_cmd=r"""java "-Xmx60g" -jar /home/kyle.lesack1/miniconda3/envs/picardtools/share/picard-2.27.5-0/picard.jar SortSam """,
#			max_records="25000"
#	conda:  "yaml/picardtools.v2.27.5.yaml"
#	threads: 8
#	resources:
#	        mem_mb=lambda _, attempt: 60000 + ((attempt - 1) * 10000),
#	        time_hms="02:00:00"
#	shell:
#	        "{params.picard_cmd} -I {input} -O {output.bamfile} -SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT --CREATE_INDEX --MAX_RECORDS_IN_RAM {params.max_records}"
