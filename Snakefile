ALL_STRAINS = ['DL238', 'ECA36', 'ECA396', 'EG4725', 'JU1400', 'JU2526', 'JU2600', 'JU310', 'MY2147', 'MY2693', 'NIC2', 'NIC526', 'QX1794', 'XZ1516']
REFERENCE = "0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa" # C. elegans reference genome
NGMLRDICT40X = {"DL238": "0.24", "DRR142768": "0.52", "ECA36": "0.19", "ECA396": "0.30", "EG4725": "0.18", "JU1400": "0.34", "JU2526": "0.33", "JU2600": "0.22", "JU310": "0.25", "MY2147":"0.21", "MY2693": "0.31", "NIC2": "0.33", "NIC526": "0.27", "QX1794": "0.32","XZ1516": "0.54"}
ALIGNMENT_DIR = ["ngmlr", "ngmlr.subsampled_40X"]

rule all:
	input:
		#expand("1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.{extension}", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS, extension=["bam","bam.bai"]),
		#expand("1_alignments/ngmlr.subsampled_40X/{strain}/{strain}_picard_sorted.bam", strain=ALL_STRAINS),
		#expand("1_alignments/ngmlr.subsampled.picard_sorted/{strain}/{strain}_picard_sorted.bam", strain=ALL_STRAINS), # use diff on these and subsampled. Remove this if same.
		expand("2_variant_calls/{alignment_dir}/sniffles/{strain}/{strain}.{extension}", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS, extension=["vcf", "vcf.gz"]),
		expand("3_jasmine/{alignment_dir}/sniffles/dup_to_ins/{strain}_dupToIns.vcf", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS),
		expand("3_jasmine/{alignment_dir}/sniffles/dup_to_ins/{strain}_dupToIns_refined.vcf", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS),
		#expand("2_variant_calls/{alignment_dir}/sniffles/{strain}/variants.sorted.vcf", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS),
		#expand("2_variant_calls/{alignment_dir}/sniffles/{strain}/filtered/passed/passed.sorted.{extension}", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS, extension = ["vcf","vcf.gz"]),
		#expand("2_variant_calls/{alignment_dir}/sniffles/{strain}/filtered/passed/merged/passed.sorted.merged.vcf.{extension}", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS, extension = ["gz", "gz.tbi"]),
		#expand("2_variant_calls/{alignment_dir}/sniffles/all_samples_merged/passed.sorted.merged.vcf.gz", alignment_dir = ALIGNMENT_DIR),
		#expand("2_variant_calls/{alignment_dir}/sniffles/all_samples_merged/sv_types/svim_{svtype}.vcf", alignment_dir = ALIGNMENT_DIR, svtype = ["deletions", "tandem_duplications", "interspersed_duplications", "inversions"])

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
rule sniffles:
	input:
	        bamfile="1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.bam",
	        bamindex="1_alignments/ngmlr/{strain}/{strain}_picard_sorted.bam.bai"
	output:
		"2_variant_calls/{alignment_dir}/sniffles/{strain}/{strain}.vcf.gz"
	params:
		minsize="100",
	conda:  "yaml/sniffles2.yaml"
	threads: 8
	resources:
		mem_mb=lambda _, attempt: 20000 + ((attempt - 1) * 10000),
		time_hms="02:00:00"
	shell:
		"sniffles --input {input.bamfile} -n -1 --snf 2_variant_calls/{wildcards.alignment_dir}/sniffles/{wildcards.strain}/{wildcards.strain}.snf --minsvlen {params.minsize} --output-rnames --reference {REFERENCE} --sample-id {wildcards.strain} -t 8 --vcf {output} "

# Call SVs with SVIM
rule sort_svim:
	input:
		"2_variant_calls/{alignment_dir}/sniffles/{strain}/variants.vcf"
	output:
		"2_variant_calls/{alignment_dir}/sniffles/{strain}/variants.sorted.vcf"
	params:
	conda:  "yaml/picardtools.v2.27.5.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:10:00"
	shell:
		"picard SortVcf I={input} O={output}"

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

rule filter_svim:
	input:
		"2_variant_calls/{alignment_dir}/sniffles/{strain}/variants.sorted.vcf"
	output:
		expand("2_variant_calls/{alignment_dir}/sniffles/{strain}/filtered/passed/passed.sorted.{extension}", extension = ["vcf", "vcf.gz"], allow_missing=True),
	params:
		SVIM_SV_TYPES = ["INV", "DEL","DUP:INT", "DUP:TANDEM","DUP"], # I didn't include INS or BND. I'm only interested in variants that span genes.
		outdir = "2_variant_calls/{alignment_dir}/sniffles/{strain}/filtered/",
		minsize = "100",
		min_qual_svim = "15",
		sample_name = "{strain}"
	conda:  "yaml/pysam.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	script:
		"scripts/1_process_results/parse_svim.py"

rule bgzip_decompress_vcf_gz:
	input:
		"2_variant_calls/{alignment_dir}/sniffles/{strain}/{strain}.vcf.gz"
	output:
		"2_variant_calls/{alignment_dir}/sniffles/{strain}/{strain}.vcf"
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""bgzip -d -k {input}"""

rule jasmine_dup_to_ins:
	input:
		expand("2_variant_calls/{alignment_dir}/sniffles/{strain}/{strain}.vcf", strain=ALL_STRAINS, allow_missing=True),
	output:
		expand("3_jasmine/{alignment_dir}/sniffles/dup_to_ins/{strain}_dupToIns.vcf", strain=ALL_STRAINS, allow_missing=True),
	params:
		reference_genome = REFERENCE
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			echo {input} | tr " " "\n" > 3_jasmine/{wildcards.alignment_dir}/sniffles/dup_to_ins/filelist.txt
			jasmine --dup_to_ins --preprocess_only file_list=3_jasmine/{wildcards.alignment_dir}/sniffles/dup_to_ins/filelist.txt out_file={output} genome_file={params.reference_genome}
			mv output/* 3_jasmine/{wildcards.alignment_dir}/sniffles/dup_to_ins/
		"""

rule iris:
	input:
		"3_jasmine/{alignment_dir}/sniffles/dup_to_ins/{strain}_dupToIns.vcf"
	output:
		"3_jasmine/{alignment_dir}/sniffles/dup_to_ins/{strain}_dupToIns_refined.vcf"
	params:
		reference_genome = REFERENCE,
		bamfile="1_alignments/{wildcards.alignment_dir}/{wildcards.strain}/{wildcards.strain}_picard_sorted.bam",
		min_ins_length="100",
		out_dir="3_jasmine/{wildcards.alignment_dir}/sniffles/dup_to_ins/"
	conda:  "yaml/jasmine.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 5000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
		iris genome_in={params.reference_genome} vcf_in={input} reads_in={params.bamfile} vcf_out={output} min_ins_length={params.min_ins_length} out_dir={params.out_dir}
		"""

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

# Truvari Collapse
# truvari collapse -i 2_variant_calls/ngmlr/sniffles/DRR142768/filtered/passed/passed.vcf.gz -o test/merged.vcf -c test/collapsed.vcf -f 0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa --sizemin 100 --sizemax 100000


# VEP http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic

# https://github.com/brentp/duphold
