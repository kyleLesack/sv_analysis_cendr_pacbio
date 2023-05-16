import os
import sys
from collections import defaultdict
from pathlib import Path
#import argparse
#parser = argparse.ArgumentParser()
#parser.add_argument("input_vcf", help="input vcf to validate")
#parser.add_argument("output_vcf", help="validated output vcf") # The original FASTQ was used to create shuffled FASTQ files
#parser.add_argument("min_support", type=int, help="minimum variant read support")
#parser.add_argument("min_len", type=int, help="minimum variant size in bp")
#args = parser.parse_args()

with open(snakemake.log[0], "w") as f:
#with open("logs/validate_is_specific/ngmlr.txt ", "w") as f:
	sys.stderr = sys.stdout = f

	# Parse the vcf file lines to get the calls for each type
	def parse_vcf(variants, min_support, min_len):
		vcf_header_lines = []
		vcf_lines = []
		bad_lines = []
		#sv_dict = defaultdict(lambda: defaultdict(list)) # Store all variant coordinates in dictionary using filters and sv types as keys

		for line in variants:
			if line[0] == "#":
				vcf_lines.append(line)

			else:
				line_split = line.split()
				info_line = line_split[7]
				info_line_split = info_line.split(";")
				is_specific = None
				read_support = None
				sv_len = None

				for x in info_line_split: # Get the variant type
					if "SUPPORT=" in x:
						read_support = x.split("=")[1]
					elif "IS_SPECIFIC=" in x:
						is_specific = x.split("=")[1]
					elif "SVLEN=" in x:
						sv_len = x.split("=")[1]
						sv_len = abs(int(sv_len))

				if is_specific is not None:
					if read_support is not None:
						if int(read_support) >= min_support:
							if sv_len is not None:
								if int(sv_len) >= min_len:
									if is_specific != 1:
										#print("Incorrect IS_SPECIFIC value for: ")
										#print(line)
										bad_lines.append(line)
										fixed_line = line.replace("IS_SPECIFIC=0", "IS_SPECIFIC=1")
										vcf_lines.append(fixed_line)
									else:
										vcf_lines.append(line)
								else:
									vcf_lines.append(line)
							else:
								print("Could not find SVLEN for:")
								print(line)
						else:
							vcf_lines.append(line)
					else:
						print("SUPPORT not found for the following line:")
						print(line)

				else:
					print("IS_SPECIFIC not found for the following line:")
					print(line)

		return (vcf_lines, bad_lines)

	# Write to disk
	def write_vcf(vcf_lines, outdir, vcf_file):
		# Check if output directories exist.
		#outdir = os.path.dirname(outdir)
		if not os.path.exists(outdir):
			print("Creating directory: " + str(outdir))
			os.makedirs(outdir)

		print("Writing to: " + vcf_file)
		with open(vcf_file, 'w', newline='\n') as f:
			f.writelines("%s\n" % l.rstrip() for l in vcf_lines)


	min_support = snakemake.params['min_support']
	min_len = snakemake.params['min_len']
	#min_support = args.min_support
	#min_len = args.min_len
	#with open(args.input_vcf) as f:
	vcf_file_out = snakemake.output[0]
	with open(snakemake.input[0]) as f:
		vcf_lines = f.readlines()
		variants_parsed = parse_vcf(vcf_lines, min_support, min_len)
		new_vcf_lines = variants_parsed[0]
		fixed_vcf_lines = variants_parsed[1]

		#vcf_file = args.output_vcf
		outdir = Path(vcf_file_out)
		outdir = outdir.parent.absolute()
		write_vcf(new_vcf_lines, outdir, vcf_file_out)
		bad_vcf_lines_file = vcf_file_out.replace(".vcf", ".bad.vcf")
		write_vcf(fixed_vcf_lines, outdir, bad_vcf_lines_file)

		#write_vcf(passed_vcf_lines, passed_dir, "passed.sorted.vcf")
