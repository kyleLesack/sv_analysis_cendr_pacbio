import os
import sys
import numpy as np
import upsetplot as upset
from matplotlib import pyplot
#from upsetplot import from_contents
from collections import defaultdict

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_vcf", help="VCF files containing SV Calls")
parser.add_argument("min_subset_size", type=int, help="Minimum size of svs in a category to be plotted for shared variants")
parser.add_argument("min_degree", type=int, help="Minimum number of strains to share svs for them to be plotted for shared variants")
parser.add_argument("--plot", action=argparse.BooleanOptionalAction, help="Create upset plots (default just prints singleton statistics)")
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
	singleton_sv_id_dict = defaultdict(list)
	singleton_gene_dict = defaultdict(list)
	shared_sv_id_dict = defaultdict(list)
	shared_gene_dict = defaultdict(list)
	total_singletons = 0
	total_shared = 0

	for line in variants:
		if line[0] != "#":
			line_split = line.split()
			sv_id = line_split[2]
			info_line = line_split[7]
			info_line_split = info_line.split(";")
			supp_vec = None
			supp_vec_ext = None
			gene_annotation = None
			if 'WBGene' in line:
				gene_annotation_start = line.find('WBGene')
				gene_annotation = line[gene_annotation_start:]
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
						if gene_annotation is not None:
							strain_gene_dict[strain].append(gene_annotation)
					if len(strains_with_sv) == 1:
						total_singletons += 1
						singleton_sv_id_dict[strain].append(sv_id)
						if gene_annotation is not None:
							singleton_gene_dict[strain].append(gene_annotation)
					else:
						total_shared += 1
						for strain in strains_with_sv:
							shared_sv_id_dict[strain].append(sv_id)
							if gene_annotation is not None:
								shared_gene_dict[strain].append(gene_annotation)
				else:
					print("Mismatch between SUPP_VEC and SUPP_EXT")
					print("SUPP_VEC: " + supp_vec)
					print("SUPP_EXT: " + supp_vec_ext)

	print("Total singletons: " + str(total_singletons))
	print("Total shared: " + str(total_shared))
	return (strain_sv_id_dict, strain_gene_dict, singleton_sv_id_dict, singleton_gene_dict, shared_sv_id_dict, shared_gene_dict)

with open(args.input_vcf) as f:
	vcf_lines = f.readlines()
	strain_sv_ids_genes = parse_vcf(vcf_lines)
	strain_sv_ids = strain_sv_ids_genes[0]
	strain_genes = strain_sv_ids_genes[1]
	singleton_svs = strain_sv_ids_genes[2]

print("Singleton counts")
for strain in singleton_svs:
	singleton_total = len(singleton_svs[strain])
	#print(strain + ": " + str(singleton_total))
	print("Total number of singletons found in " + strain + ": " + str(singleton_total))

if args.plot:
	strain_sv_ids_df = upset.from_contents(strain_sv_ids)
	shared_svs_cardinality = upset.query(strain_sv_ids_df, min_degree=2, sort_by="cardinality")
	shared_svs_cardinality_grouped = shared_svs_cardinality.data.groupby(['DL238', 'ECA36', 'ECA396', 'EG4725', 'JU1400', 'JU2526', 'JU2600', 'JU310', 'MY2147', 'MY2693', 'NIC2', 'NIC526', 'QX1794', 'XZ1516']).size().sort_values(ascending=False)
	shared_svs_cardinality_grouped_n20 = shared_svs_cardinality_grouped[:20]
	upset.UpSet(shared_svs_cardinality_grouped, show_counts=True, sort_by='cardinality', min_subset_size = args.min_subset_size, totals_plot_elements = 20).plot()

	shared_svs_degree = upset.query(strain_sv_ids_df, min_degree=2, sort_by="-degree")
	shared_svs_degree_grouped = shared_svs_degree.data.groupby(['DL238', 'ECA36', 'ECA396', 'EG4725', 'JU1400', 'JU2526', 'JU2600', 'JU310', 'MY2147', 'MY2693', 'NIC2', 'NIC526', 'QX1794', 'XZ1516']).size().sort_values(ascending=False)
	upset.UpSet(shared_svs_degree_grouped, show_counts=True, sort_by='-degree', min_degree=int(args.min_degree)).plot()

	pyplot.show()
