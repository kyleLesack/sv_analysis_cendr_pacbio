from collections import defaultdict
from pathlib import Path
from builtins import any as b_any

high_impact_genes = set(line.strip() for line in open(snakemake.input[0]))

gene_dict = {}
with open(snakemake.input[1], 'r') as f:
	for line in f:
		line_split = line.strip().split(" ")
		entrez_gene_name = line_split[0]
		wbgene_name = line_split[1].strip()
		
		if b_any(wbgene_name in x for x in high_impact_genes):
			gene_dict[wbgene_name]= entrez_gene_name
		else:
			print("Excluding " + wbgene_name + " due to it not having a high impact variant")

sv_names = set()
with open(snakemake.input[2], 'r') as f:
	for line in f:
		for wbgene_name in gene_dict.keys():
			if wbgene_name in line:
				line_split = line.strip().split("\t")
				sv_name = line_split[7]
				sv_names.add(sv_name)

sv_genes = defaultdict(list)			
sv_strains = defaultdict(set)			
with open(snakemake.input[2], 'r') as f:
	for line in f:
		for svname in sv_names:
			if svname in line:
				line_split = line.strip().split("\t")
				gene_spanned = line_split[3]
				gene_spanned_svname=gene_spanned + "," + svname
				if gene_spanned_svname in high_impact_genes:				
					sv_genes[svname].append(gene_spanned)
					strains_with_sv = line_split[8].strip()
					sv_strains[svname].add(strains_with_sv)
				else:
					print("Excluding " + gene_spanned_svname + " due to it not having a high impact variant")					

summary_lines = []
for sv_name in sv_genes.keys():
	genes = ','.join(sv_genes[sv_name])
	strains = ','.join(sv_strains[sv_name])
	sv_summary_line = sv_name + "\t" + genes + "\t" + strains
	summary_lines.append(sv_summary_line)

summary_lines.sort()
summary_lines.insert(0, "SVName\tGenes\tStrains")

with open(snakemake.output[0], 'w') as f:
	f.writelines(f'{s}\n' for s in summary_lines)
