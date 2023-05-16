import os
import sys
from collections import defaultdict

with open(snakemake.log[0], "w") as f:
	sys.stderr = sys.stdout = f

# Parse the vcf file lines to get the calls for each type
	def parse_vcf(variants):
		vcf_header_lines = []
		#vcf_lines = []
		#sv_dict = defaultdict(lambda: defaultdict(list)) # Store all variant coordinates in dictionary using filters and sv types as keys
		sv_dict = defaultdict(list) # Store all variant coordinates in dictionary using filters and sv types as keys

		for line in variants:
			if line[0] == "#":
				vcf_header_lines.append(line)

			else:
				line_split = line.split()
				info_line = line_split[7]
				info_line_split = info_line.split(";")
				variant_type = None

				for x in info_line_split: # Get the variant type
					if "SVTYPE=" in x:
						variant_type = x.split("=")[1]
				if variant_type is not None:
					sv_dict[variant_type].append(line)
				else:
					print("Variant type not found for the following line:")
					print(line)

		return (vcf_header_lines, sv_dict)

	# Write svs that were excluded to disk
	def write_vcf(vcf_lines, outdir, filename):
		# Check if output directories exist.
		#outdir = os.path.dirname(outdir)
		if not os.path.exists(outdir):
			print("Creating directory: " + outdir)
			os.makedirs(outdir)

		output_file = outdir + "/" + filename
		print("Writing to: " + output_file)
		with open(output_file, 'w', newline='\n') as f:
			f.writelines("%s\n" % l.rstrip() for l in vcf_lines)

	#with open("3_jasmine/ngmlr/sniffles/8_final/merged_final_genotyped.vcf") as f:

	# Process unfiltered Jasmine results
	#with open(snakemake.input[0]) as f:
	with open(snakemake.input['unfiltered']) as f:
		vcf_lines = f.readlines()
		variants_parsed = parse_vcf(vcf_lines)
		vcf_header = variants_parsed[0]
		variant_dict = variants_parsed[1]
		#outdir = os.path.dirname(snakemake.params['out_dir_unfiltered'])
		outdir = snakemake.params['out_dir_unfiltered']
		#outdir = "3_jasmine/ngmlr/sniffles/8_final/sv_types/"

		for sv_type in variant_dict.keys():
				#variant_dir = outdir
				vcf_lines = vcf_header + variant_dict[sv_type]
				vcf_file = sv_type + ".vcf"
				write_vcf(vcf_lines, outdir, vcf_file)
		#write_vcf(passed_vcf_lines, passed_dir, "passed.sorted.vcf")

	# Process filtered Jasmine results
	with open(snakemake.input['filtered']) as f:
		vcf_lines = f.readlines()
		variants_parsed = parse_vcf(vcf_lines)
		vcf_header = variants_parsed[0]
		variant_dict = variants_parsed[1]
		outdir = snakemake.params['out_dir_filtered']

		for sv_type in variant_dict.keys():
				vcf_lines = vcf_header + variant_dict[sv_type]
				vcf_file = sv_type + ".vcf"
				write_vcf(vcf_lines, outdir, vcf_file)
