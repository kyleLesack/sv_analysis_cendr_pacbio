import os
import sys
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_vcf", help="VCF files containing VEP annotations")
args = parser.parse_args()


# Parse the vcf file lines to get the calls for each type
def parse_vcf(variants):
	gene_set = set()
	gene_count = 0
	genes_affected_by_multiple_svs = 0

	for line in variants:

		if line[0] != "#":
			if "WBGene" in line:
				gene_count += line.count('WBGene')
			if line.count('WBGene') == 1:
				gene_annotation_start = line.find('WBGene')
				gene_annotation = line[gene_annotation_start:]
				#print(gene_annotation.rstrip())
				gene_annotation_end = gene_annotation.find('|')
				gene_annotation = gene_annotation[0:gene_annotation_end]
				if gene_annotation in gene_set:
					genes_affected_by_multiple_svs += 1
				gene_set.add(gene_annotation)

				#print(gene_annotation)
				#print(line[0:s.find('.')])
			else:
				print("Found multiple gene annotations")

			#line_split = line.split()
			#info_line = line_split[7]
			#info_line_split = info_line.split(";")
			#variant_type = None

			#for x in info_line_split: # Go though info field and find the VEP consequence annotations
			#	if "CSQ=" in x:
			#		csq_annotations = x.split("=")[1]
			#		for annotation in csq_annotations:
			#			print(annotation)
			#			if "WBGene" in annotation:
			#				print(annotation)

	print("# of WormBase genes in VEP results: " + str(gene_count))
	print("# of unique WormBase genes in VEP results: " + str(len(gene_set)))
	#print(genes_affected_by_multiple_svs)
	#print(len(gene_set) + genes_affected_by_multiple_svs)
	return (gene_set)

with open(args.input_vcf) as f:
	vcf_lines = f.readlines()
	gene_set = parse_vcf(vcf_lines)

write_genes(gene_set)

	#write_vcf(passed_vcf_lines, passed_dir, "passed.sorted.vcf")
