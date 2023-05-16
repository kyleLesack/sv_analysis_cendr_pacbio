import os
import sys
import numpy as np
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_vcf", help="VCF files containing VEP annotations")
args = parser.parse_args()

STRAINS = np.array(["DL238", "ECA36", "ECA396", "EG4725", "JU1400", "JU2526", "JU2600", "JU310", "MY2147", "MY2693", "NIC2", "NIC526", "QX1794", "XZ1516"])

def get_strains_with_variant(supp_vec):
	bool_list = [x=="1" for x in supp_vec]
	strains_with_sv = STRAINS[bool_list]
	return(strains_with_sv)

# Parse the vcf file lines to get the calls for each type
def parse_vcf(variants):
	strain_sv_id_dict = defaultdict(list) # Store all sv calls in dictionary using the strains as keys
	strain_gene_dict = defaultdict(list) # Store all genes impacted by sv calls in dictionary using the strains as keys

	for line in variants:
		if line[0] != "#":
			line_split = line.split()
			sv_id = line_split[2]
			info_line = line_split[7]
			info_line_split = info_line.split(";")
			supp_vec = None
			supp_vec_ext = None

			gene_annotation_start = line.find('WBGene')
			gene_annotation = line[gene_annotation_start:]
			#print(gene_annotation.rstrip())
			gene_annotation_end = gene_annotation.find('|')
			gene_annotation = gene_annotation[0:gene_annotation_end]

			for x in info_line_split: # Go though info field and find the genotype information
				if "SUPP_VEC=" in x:
					supp_vec = x.split("=")[1]
				elif "SUPP_VEC_EXT=" in x:
					supp_vec_ext = x.split("=")[1]
			if supp_vec is not None and supp_vec_ext is not None:
				if supp_vec == supp_vec_ext:
					strains_with_sv = get_strains_with_variant(supp_vec)
					for strain in strains_with_sv:
						strain_sv_id_dict[strain].append(sv_id)
						strain_gene_dict[strain].append(gene_annotation)

				else:
					print("Mismatch between SUPP_VEC and SUPP_EXT")
					print("SUPP_VEC: " + supp_vec)
					print("SUPP_EXT: " + supp_vec_ext)

	return (strain_sv_id_dict, strain_gene_dict)

with open(args.input_vcf) as f:
	vcf_lines = f.readlines()
	strain_sv_ids_genes = parse_vcf(vcf_lines)
	strain_sv_ids = strain_sv_ids_genes[0]
	strain_genes = strain_sv_ids_genes[1]

for strain in strain_sv_ids:
	sv_total = len(strain_sv_ids[strain])
	print("Total number of variants found in " + strain + ": " + str(sv_total))
