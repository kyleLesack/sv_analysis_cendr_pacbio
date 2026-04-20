from collections import defaultdict
from pathlib import Path
from statistics import mean, median, mode

sv_genes = defaultdict(list)
with open(snakemake.input[0], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t")
			gene_name = line_split[3]
			sv_name = line_split[7]
			sv_genes[sv_name].append(gene_name)

gene_count = []
overlapped_genes = []
most_overlap_count = 0
most_overlap_sv = None
most_overlap_genes = None
for sv in sv_genes.keys():
	gene_list = sv_genes[sv]
	gene_count.append(len(gene_list))
	for gene in gene_list:
		overlapped_genes.append(gene)
	if len(gene_list) > most_overlap_count:
		most_overlap_count = len(gene_list)
		most_overlap_sv = sv
		most_overlap_genes = gene_list

# Store summary statistics in list that will be written to disk
summary_lines = []
summary_lines.append("Highest overlap count: " + str(most_overlap_count))
summary_lines.append("Highest overlapping sv:" + most_overlap_sv)
summary_lines.append(most_overlap_genes)
summary_lines.append("Summary statistics for the number of genes\n")
summary_lines.append("Genes overlapped by SVs: " + str(len(set(overlapped_genes))))
summary_lines.append("Mean number of genes overlapped: " + str((round(mean(gene_count),2))))
summary_lines.append("Mode of genes overlapped: " + str(((mode(gene_count)))))
summary_lines.append("Minimum number of genes overlapped: " + str(min(gene_count)))
summary_lines.append("Maximum number of genes overlapped: " + str(max(gene_count)))

with open(snakemake.output[0], 'w') as f:
	f.writelines(f'{s}\n' for s in summary_lines)