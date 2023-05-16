ALL_STRAINS = ['DL238', 'ECA36', 'ECA396', 'EG4725', 'JU1400', 'JU2526', 'JU2600', 'JU310', 'MY2147', 'MY2693', 'NIC2', 'NIC526', 'QX1794', 'XZ1516']
REFERENCE = "0_input/reference/c_elegans.PRJNA13758.WS263.genomic.fa" # C. elegans reference genome
#ALIGNMENT_DIR = ["ngmlr", "ngmlr.subsampled_40X"]
ALIGNMENT_DIR = ["ngmlr"]
MIN_SUPPORT_JASMINE = {"DL238": "42", "ECA36": "51", "ECA396": "33", "EG4725": "54", "JU1400": "29", "JU2526": "30", "JU2600": "44", "JU310": "40", "MY2147": "46", "MY2693": "32", "NIC2": "30", "NIC526": "37", "QX1794": "31", "XZ1516": "18"} # 25% of the sequencing depth for each sample bam file

rule all:
	input:
		#expand("1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.{extension}", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS, extension=["bam","bam.bai"]),
		#expand("1_alignments/ngmlr.subsampled_40X/{strain}/{strain}_picard_sorted.bam", strain=ALL_STRAINS),
		#expand("1_alignments/ngmlr.subsampled.picard_sorted/{strain}/{strain}_picard_sorted.bam", strain=ALL_STRAINS), # use diff on these and subsampled. Remove this if same.
		#expand("2_variant_calls/{alignment_dir}/sniffles/{strain}/{strain}.{extension}", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS, extension=["vcf", "vcf.gz"]),
		expand("3_jasmine/{alignment_dir}/sniffles/1_dup_to_ins/{strain}_dupToIns.vcf", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS),
		expand("3_jasmine/{alignment_dir}/sniffles/2_iris_refined/{strain}_dupToIns_refined.vcf", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS),
		expand("3_jasmine/{alignment_dir}/sniffles/3_normalized/{strain}_dupToIns_refined_normalizeTypes.vcf", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS),
		expand("3_jasmine/{alignment_dir}/sniffles/4_specific/{strain}_dupToIns_refined_normalizeTypes_markedSpec.vcf", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS),
		expand("3_jasmine/{alignment_dir}/sniffles/5_no_dupes/{strain}_dupToIns_refined_normalizeTypes_markedSpec_no_dupes.vcf", alignment_dir = ALIGNMENT_DIR, strain=ALL_STRAINS),
		expand("3_jasmine/{alignment_dir}/sniffles/6_merged/merged.vcf", alignment_dir = ALIGNMENT_DIR),
		expand("3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/merged_ins_to_dup.vcf", alignment_dir = ALIGNMENT_DIR),
		expand("3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/sv_types/{variant_type}.vcf", alignment_dir = ALIGNMENT_DIR, variant_type = ["DEL", "DUP", "INV"]),
		expand("3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/validated_is_specific/{filename}", alignment_dir = ALIGNMENT_DIR, filename = ["merged_ins_to_dup.vcf", "merged_ins_to_dup_genotyped.vcf"]),
		expand("3_jasmine/{alignment_dir}/sniffles/8_final/merged_final.vcf", alignment_dir = ALIGNMENT_DIR),
		expand("3_jasmine/{alignment_dir}/sniffles/8_final/sv_types/{variant_type}.vcf", alignment_dir = ALIGNMENT_DIR, variant_type = ["DEL", "DUP", "INV"]),
		expand("4_vep/scripts/{alignment_dir}/{input_dir}/{variant_type}-vep_per_gene-pick-no_intergenic.sh", alignment_dir = ALIGNMENT_DIR, input_dir = ["7_ins_to_dup", "8_final"], variant_type = ["DEL", "DUP", "INV"]),
		expand("4_vep/scripts/{alignment_dir}/{input_dir}/annotated_high_impact/{variant_type}-vep_per_gene-pick-no_intergenic-annotated_high_impact.sh", alignment_dir = ALIGNMENT_DIR, input_dir = ["7_ins_to_dup", "8_final"], variant_type = ["DEL", "DUP", "INV"]),
		#expand("5_wormcat/scripts/{alignment_dir}/{input_dir}/{variant_type}/{variant_type}-run_wormcat_on_vep_results.sh", alignment_dir = ALIGNMENT_DIR, input_dir = ["7_ins_to_dup", "8_final"], variant_type = ["DEL", "DUP", "INV"])
		expand("5_wormcat/output/{alignment_dir}/{input_dir}/{variant_type}/{variant_type}-vep_genes_high_impact_wormcat.R", alignment_dir = ["ngmlr"], input_dir = ["7_ins_to_dup", "8_final"], variant_type = ["DEL", "DUP", "INV"])

rule picard_sort:
	input:
			"1_alignments/{alignment_dir}/{strain}/{strain}_sorted.bam"
	output:
	        bamfile="1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.bam",
	        bamindex="1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.bai"
	params:
	        picard_cmd=r"""java "-Xmx60g" -jar /home/kyle.lesack1/miniconda3/envs/picardtools/share/picard-2.27.5-0/picard.jar SortSam """,
			max_records="25000",
	conda:  "yaml/picardtools.v2.27.5.yaml"
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

# Call SVs with Sniffles
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
		"sniffles --input {input.bamfile} --snf 2_variant_calls/{wildcards.alignment_dir}/sniffles/{wildcards.strain}/{wildcards.strain}.snf --minsvlen {params.minsize} --output-rnames --reference {REFERENCE} --sample-id {wildcards.strain} -t 8 --vcf {output} "

# Decompress gz files
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

# Convert duplications in Jasmine in order to refine them using Iris
rule jasmine_dup_to_ins:
	input:
		expand("2_variant_calls/{alignment_dir}/sniffles/{strain}/{strain}.vcf", strain=ALL_STRAINS, allow_missing=True),
	output:
		expand("3_jasmine/{alignment_dir}/sniffles/1_dup_to_ins/{strain}_dupToIns.vcf", strain=ALL_STRAINS, allow_missing=True),
	params:
		reference_genome = REFERENCE,
		out_dir="3_jasmine/{alignment_dir}/sniffles/1_dup_to_ins/",
		out_file="3_jasmine/{alignment_dir}/sniffles/1_dup_to_ins/dummy.vcf" # I don't think this is used when Jasmine is run with the option --preprocess_only
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			echo {input} | tr " " "\n" > 3_jasmine/{wildcards.alignment_dir}/sniffles/1_dup_to_ins/sniffles_vcf_files.txt
			jasmine --dup_to_ins --preprocess_only file_list=3_jasmine/{wildcards.alignment_dir}/sniffles/1_dup_to_ins/sniffles_vcf_files.txt out_file={output} genome_file={params.reference_genome} out_dir={params.out_dir}
		"""

# Refine insertions with Iris
rule iris:
	input:
		"3_jasmine/{alignment_dir}/sniffles/1_dup_to_ins/{strain}_dupToIns.vcf"
	output:
		"3_jasmine/{alignment_dir}/sniffles/2_iris_refined/{strain}_dupToIns_refined.vcf"
	params:
		reference_genome = REFERENCE,
		bamfile="1_alignments/{alignment_dir}/{strain}/{strain}_picard_sorted.bam",
		min_ins_length="100",
		out_dir="3_jasmine/{alignment_dir}/sniffles/2_iris_refined/iris_intermediate_files/{strain}"
	conda:  "yaml/iris.yaml"
	threads: 4
	resources:
		mem_mb=lambda _, attempt: 20000 + ((attempt - 1) * 10000),
		time_hms="05:00:00"
	shell:
		"""
			mkdir -p {params.out_dir}
			iris genome_in={params.reference_genome} vcf_in={input} reads_in={params.bamfile} vcf_out={output} min_ins_length={params.min_ins_length} out_dir={params.out_dir}
		"""
# Documented as "run type normalization before merging". Seems to convert BND to TRA
rule jasmine_normalize:
	input:
		expand("3_jasmine/{alignment_dir}/sniffles/2_iris_refined/{strain}_dupToIns_refined.vcf", strain=ALL_STRAINS, allow_missing=True),
	output:
		expand("3_jasmine/{alignment_dir}/sniffles/3_normalized/{strain}_dupToIns_refined_normalizeTypes.vcf", strain=ALL_STRAINS, allow_missing=True),
	params:
		reference_genome = REFERENCE,
		out_dir="3_jasmine/{alignment_dir}/sniffles/3_normalized/",
		out_file="3_jasmine/{alignment_dir}/sniffles/3_normalized/dummy.vcf"
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			echo {input} | tr " " "\n" > 3_jasmine/{wildcards.alignment_dir}/sniffles/3_normalized/sniffles_vcf_files_refined.txt
			jasmine --pre_normalize  --preprocess_only file_list=3_jasmine/{wildcards.alignment_dir}/sniffles/3_normalized/sniffles_vcf_files_refined.txt out_file={params.out_file} genome_file={params.reference_genome} out_dir={params.out_dir}
		"""

# Documented as "mark calls in the original VCF files that have enough support to called specific". On GitHub pipeline it is described as "Mark high-confidence callset (high-specificity callset) in each sample"
rule jasmine_mark_specific:
	input:
		"3_jasmine/{alignment_dir}/sniffles/3_normalized/{strain}_dupToIns_refined_normalizeTypes.vcf"
	output:
		"3_jasmine/{alignment_dir}/sniffles/4_specific/{strain}_dupToIns_refined_normalizeTypes_markedSpec.vcf"
	wildcard_constraints:
		alignment_dir = "ngmlr"
	params:
		reference_genome = REFERENCE,
		out_dir="3_jasmine/{alignment_dir}/sniffles/4_specific/",
		spec_len="100",
		spec_reads=lambda wcs: MIN_SUPPORT_JASMINE[wcs.strain],
		file_list="3_jasmine/{alignment_dir}/sniffles/4_specific/{strain}_sniffles_vcf_file_refined.txt",
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			echo {input} | tr " " "\n" > {params.file_list}
			jasmine --mark_specific spec_reads={params.spec_reads} spec_len={params.spec_len} --preprocess_only file_list={params.file_list} out_file={output} genome_file={params.reference_genome} out_dir={params.out_dir}
		"""

rule jasmine_mark_specific_40x:
	input:
		"3_jasmine/{alignment_dir}/sniffles/3_normalized/{strain}_dupToIns_refined_normalizeTypes.vcf"
	output:
		"3_jasmine/{alignment_dir}/sniffles/4_specific/{strain}_dupToIns_refined_normalizeTypes_markedSpec.vcf"
	wildcard_constraints:
		alignment_dir = "ngmlr.subsampled_40X"
	params:
		reference_genome = REFERENCE,
		out_dir="3_jasmine/{alignment_dir}/sniffles/4_specific/",
		spec_len="100",
		spec_reads="10",
		file_list="3_jasmine/{alignment_dir}/sniffles/4_specific/{strain}_sniffles_vcf_file_refined.txt",
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			echo {input} | tr " " "\n" > {params.file_list}
			jasmine --mark_specific spec_len={params.spec_len} spec_reads=10 --preprocess_only file_list={params.file_list} out_file={output} genome_file={params.reference_genome} out_dir={params.out_dir}
		"""

# Remove duplicate calls in each sample
rule jasmine_remove_duplicates:
	input:
		"3_jasmine/{alignment_dir}/sniffles/4_specific/{strain}_dupToIns_refined_normalizeTypes_markedSpec.vcf"
	output:
		"3_jasmine/{alignment_dir}/sniffles/5_no_dupes/{strain}_dupToIns_refined_normalizeTypes_markedSpec_no_dupes.vcf"
	params:
		reference_genome = REFERENCE,
		out_dir="3_jasmine/{alignment_dir}/sniffles/5_no_dupes/",
		file_list="3_jasmine/{alignment_dir}/sniffles/5_no_dupes/{strain}_sniffles_vcf_files_specific.txt",
		max_dist="200",
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			echo {input} > {params.file_list}
			jasmine max_dist={params.max_dist} --allow_intrasample --nonlinear_dist file_list={params.file_list} out_file={output} genome_file={params.reference_genome} out_dir={params.out_dir}
		"""

# Merge all samples into single vcf
rule jasmine_merge_samples:
	input:
		expand("3_jasmine/{alignment_dir}/sniffles/5_no_dupes/{strain}_dupToIns_refined_normalizeTypes_markedSpec_no_dupes.vcf", strain = ALL_STRAINS, allow_missing = True)
	output:
		"3_jasmine/{alignment_dir}/sniffles/6_merged/merged.vcf"
	params:
		reference_genome = REFERENCE,
		out_dir="3_jasmine/{alignment_dir}/sniffles/6_merged/",
		file_list="3_jasmine/{alignment_dir}/sniffles/6_merged/sniffles_vcf_file_specific_no_dupes.txt",
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			echo {input} | tr " " "\n" > {params.file_list}
			jasmine --output_genotypes file_list={params.file_list} out_file={output} genome_file={params.reference_genome} out_dir={params.out_dir}
		"""

# Convert insertions back to duplications
rule jasmine_ins_to_dup:
	input:
		"3_jasmine/{alignment_dir}/sniffles/6_merged/merged.vcf"
	output:
		genotyped="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/merged_ins_to_dup_genotyped.vcf",
		not_genotyped="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/merged_ins_to_dup.vcf"
	params:
		out_dir="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/",
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			cp {input} {output.genotyped}
			cp {input} {output.not_genotyped}
			cp {input} {output.genotyped}.bak
			cp {input} {output.not_genotyped}.bak
			jasmine --dup_to_ins --postprocess_only out_file={output.genotyped} out_dir={params.out_dir}
			jasmine --dup_to_ins --postprocess_only out_file={output.not_genotyped} out_dir={params.out_dir}
		"""

# Check if Jasmine added the correct values for IS_SPECIFIC for the subsampled bams
rule validate_is_specific_40x:
	input:
		"3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/{filename}"
	output:
		"3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/validated_is_specific/{filename}"
	wildcard_constraints:
		alignment_dir = "ngmlr.subsampled_40X"
	log:
		"logs/validate_is_specific/{alignment_dir}_{filename}.txt"
	params:
		min_support = 10,
		min_len = 100
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	script:
		"scripts/1_process_results/check_jasmin_specific_calls.py"

# Check if Jasmine added the correct values for IS_SPECIFIC for the full depth bams
rule validate_is_specific:
	input:
		"3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/{filename}"
	output:
		"3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/validated_is_specific/{filename}"
	wildcard_constraints:
		alignment_dir = "ngmlr"
	log:
		"logs/validate_is_specific/{alignment_dir}_{filename}.txt"
	params:
		min_support = 18, # lowest value in the MIN_SUPPORT_JASMINE dict
		min_len = 100
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	script:
		"scripts/1_process_results/check_jasmin_specific_calls.py"

# Finalize merged dataset. Remove low confidence calls.
rule jasmine_finalize_40x:
	input:
		genotyped="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/validated_is_specific/merged_ins_to_dup_genotyped.vcf",
		not_genotyped="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/validated_is_specific/merged_ins_to_dup.vcf"
	output:
		genotyped="3_jasmine/{alignment_dir}/sniffles/8_final/merged_final_genotyped.vcf",
		not_genotyped="3_jasmine/{alignment_dir}/sniffles/8_final/merged_final.vcf"
	wildcard_constraints:
		alignment_dir = "ngmlr.subsampled_40X"
	params:
		out_dir="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/",
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			cat {input.genotyped} | grep -v 'IMPRECISE;' | grep -v 'IS_SPECIFIC=0' > {output.genotyped}
			cat {input.not_genotyped} | grep -v 'IMPRECISE;' | grep -v 'IS_SPECIFIC=0' > {output.not_genotyped}
		"""

# Finalize merged dataset. Remove low confidence calls.
rule jasmine_finalize:
	input:
		genotyped="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/merged_ins_to_dup_genotyped.vcf",
		not_genotyped="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/merged_ins_to_dup.vcf"
	output:
		genotyped="3_jasmine/{alignment_dir}/sniffles/8_final/merged_final_genotyped.vcf",
		not_genotyped="3_jasmine/{alignment_dir}/sniffles/8_final/merged_final.vcf"
	wildcard_constraints:
		alignment_dir = "ngmlr"
	params:
		out_dir="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/",
	conda:  "yaml/jasmine.yaml"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	shell:
		"""
			cat {input.genotyped} | grep -v 'IMPRECISE;' | grep -v 'IS_SPECIFIC=0' > {output.genotyped}
			cat {input.not_genotyped} | grep -v 'IMPRECISE;' | grep -v 'IS_SPECIFIC=0' > {output.not_genotyped}
		"""

# Create different files for the Jasmine results for each variant type
rule extract_svs_by_type_40x:
	input:
		unfiltered="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/validated_is_specific/merged_ins_to_dup.vcf",
		filtered="3_jasmine/{alignment_dir}/sniffles/8_final/merged_final.vcf"
	output:
		unfiltered=expand("3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/sv_types/{variant_type}.vcf", variant_type = ["DEL", "DUP", "INV"], allow_missing = True),
		filtered=expand("3_jasmine/{alignment_dir}/sniffles/8_final/sv_types/{variant_type}.vcf", variant_type = ["DEL", "DUP", "INV"], allow_missing = True)
	wildcard_constraints:
		alignment_dir = "ngmlr.subsampled_40X"
	log:
		"logs/extract_svs_by_type/{alignment_dir}.txt"
	params:
		out_dir_unfiltered = "3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/sv_types/",
		out_dir_filtered = "3_jasmine/{alignment_dir}/sniffles/8_final/sv_types/",
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	script:
		"scripts/1_process_results/extract_svs_by_type.py"

# Create different files for the Jasmine results for each variant type
rule extract_svs_by_type:
	input:
		unfiltered="3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/merged_ins_to_dup.vcf",
		filtered="3_jasmine/{alignment_dir}/sniffles/8_final/merged_final.vcf"
	output:
		unfiltered=expand("3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/sv_types/{variant_type}.vcf", variant_type = ["DEL", "DUP", "INV"], allow_missing = True),
		filtered=expand("3_jasmine/{alignment_dir}/sniffles/8_final/sv_types/{variant_type}.vcf", variant_type = ["DEL", "DUP", "INV"], allow_missing = True),
	wildcard_constraints:
		alignment_dir = "ngmlr"
	log:
		"logs/extract_svs_by_type/{alignment_dir}.txt"
	params:
		out_dir_unfiltered = "3_jasmine/{alignment_dir}/sniffles/7_ins_to_dup/sv_types/",
		out_dir_filtered = "3_jasmine/{alignment_dir}/sniffles/8_final/sv_types/",
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:05:00"
	script:
		"scripts/1_process_results/extract_svs_by_type.py"


# Create scripts for annotating variants with VEP, as Docker doesn't play nice with Snakemake
rule vep_annotate:
	input:
		"3_jasmine/{alignment_dir}/sniffles/{input_dir}/sv_types/{variant_type}.vcf"
	output:
		"4_vep/scripts/{alignment_dir}/{input_dir}/{variant_type}-vep_per_gene-pick-no_intergenic.sh"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:20:00"
	shell:
		"""
			echo mkdir -p ./4_vep/input/{wildcards.alignment_dir}/{wildcards.input_dir}/ > {output}
			echo mkdir -p ./4_vep/output/{wildcards.alignment_dir}/{wildcards.input_dir}/ >> {output}
			echo cp {input} ./4_vep/input/{wildcards.alignment_dir}/{wildcards.input_dir}/ >> {output}
			echo chmod -R a+rwx ./4_vep >> {output}
			echo docker run -v /bulk/worm_lab/mrkyle/sv_analysis_cendr_pacbio/4_vep:/data ensemblorg/ensembl-vep \
			  vep --cache --offline --format vcf --vcf --force_overwrite \
			      --input_file input/{wildcards.alignment_dir}/{wildcards.input_dir}/{wildcards.variant_type}.vcf \
			      --output_file output/{wildcards.alignment_dir}/{wildcards.input_dir}/{wildcards.variant_type}-vep_per_gene-pick-no_intergenic.vcf \
			      --species caenorhabditis_elegans \
			      --per_gene \
			      --pick \
			      --no_intergenic >> {output}
		"""

# Create scripts to process VEP results by removing VCF lines without annotations
rule filter_vep_annotated_high_impact:
	input:
		"4_vep/scripts/{alignment_dir}/{input_dir}/{variant_type}-vep_per_gene-pick-no_intergenic.sh"
	output:
		"4_vep/scripts/{alignment_dir}/{input_dir}/annotated_high_impact/{variant_type}-vep_per_gene-pick-no_intergenic-annotated_high_impact.sh"
	params:
		tempdir = "4_vep/output/{alignment_dir}/{input_dir}/{variant_type}_temp",
		docker_tempdir = "output/{alignment_dir}/{input_dir}/{variant_type}_temp",
		grep_cmd = r"""grep "#" """,
		grep_invert_cmd = r"""grep -v "#" """,
		impact_setting = r"""filter "IMPACT is HIGH" """
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:20:00"
	shell:
		"""
			echo mkdir -p {params.tempdir} > {output}
			echo {params.grep_cmd:q} 4_vep/output/{wildcards.alignment_dir}/{wildcards.input_dir}/{wildcards.variant_type}-vep_per_gene-pick-no_intergenic.vcf \> {params.tempdir}/{wildcards.variant_type}-vep_vep_per_gene-pick-no_intergenic_annotated_only.vcf >> {output}
			echo {params.grep_invert_cmd:q} 4_vep/output/{wildcards.alignment_dir}/{wildcards.input_dir}/{wildcards.variant_type}-vep_per_gene-pick-no_intergenic.vcf \| grep "CSQ" \>\> {params.tempdir}/{wildcards.variant_type}-vep_vep_per_gene-pick-no_intergenic_annotated_only.vcf >> {output}
			echo mkdir -p ./4_vep/output/{wildcards.alignment_dir}/{wildcards.input_dir}/annotated_high_impact/ >> {output}
			echo chmod -R a+rwx ./4_vep >> {output}
			echo docker run -v /bulk/worm_lab/mrkyle/sv_analysis_cendr_pacbio/4_vep:/data ensemblorg/ensembl-vep \
			  filter_vep --format vcf --vcf --force_overwrite \
			      --input_file {params.docker_tempdir}/{wildcards.variant_type}-vep_vep_per_gene-pick-no_intergenic_annotated_only.vcf \
			      --output_file output/{wildcards.alignment_dir}/{wildcards.input_dir}/annotated_high_impact/{wildcards.variant_type}-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf \
			      --vcf_info_field CSQ \
			      --{params.impact_setting:q} >> {output}
  		  	echo rm -rf {params.tempdir} >> {output}
		"""

# Create scripts to run WormCat on VEP results
rule wormcat:
	input:
		"4_vep/scripts/{alignment_dir}/{input_dir}/annotated_high_impact/{variant_type}-vep_per_gene-pick-no_intergenic-annotated_high_impact.sh"
	output:
		"5_wormcat/output/{alignment_dir}/{input_dir}/{variant_type}/{variant_type}-vep_genes_high_impact_wormcat.R"
	params:
		vep_results = "4_vep/output/{alignment_dir}/{input_dir}/annotated_high_impact/{variant_type}-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf",
		vep_genes_high_impact = "5_wormcat\/output\/{alignment_dir}\/{input_dir}\/{variant_type}\/{variant_type}-vep_genes_high_impact.txt",
		outdir = "5_wormcat\/output\/{alignment_dir}\/{input_dir}\/{variant_type}\/wormcat_results"
		#wormcat_results= "5_wormcat/output/{alignment_dir}/{input_dir}/{variant_type}/{variant_type}-wormcat_results.txt"
	threads: 1
	resources:
		mem_mb=lambda _, attempt: 1000 + ((attempt - 1) * 10000),
		time_hms="00:20:00"
	shell:
		"""
			python3 scripts/1_process_results/extract_genes_from_vep.py {params.vep_results} {params.vep_genes_high_impact}
			sed 's/VEP_GENES_FILE/{params.vep_genes_high_impact}/g' scripts/2_wormcat/wormcat_r_cmds.R | sed 's/VARIANT_TYPE/{wildcards.variant_type}/g' | sed 's/OUTDIR/{params.outdir}/g' > {output}
			Rscript {output}
		"""


# WormCat
# Don't use the WormCat results that are corrected for multiple comparisons using Bonferonni. Use a family-wise error rate approach.
# https://cran.r-project.org/web/packages/cherry/index.html
# See Multiple Testing of Gene Sets from Gene Ontology: Possibilities and Pitfalls

# https://github.com/brentp/duphold
