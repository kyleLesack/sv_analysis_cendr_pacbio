import os
import sys
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_vcf", help="VCF file to convert to bed")
args = parser.parse_args()


# Parse the vcf file lines to get the calls for each type
def parse_vcf(variants):

	for line in variants:
		end_coord = None
		if line[0] != "#":
			line_split = line.split()
			info_line = line_split[7]
			info_line_split = info_line.split(";")
			chromosome = line_split[0]
			start_coord = line_split[1]
			for x in info_line_split: # Go though info field and find the genotype information
				if "END=" in x:
					x_split = x.split("=")
					if x_split[0] == "END":
						end_coord = x.split("=")[1]

			if end_coord is not None:
				bed_line = chromosome + "\t" + start_coord + "\t" + end_coord
				print(bed_line)
			else:
				print("Did not find end coordinate for: " + line)
				input("")



with open(args.input_vcf) as f:
	vcf_lines = f.readlines()
	parse_vcf(vcf_lines)
