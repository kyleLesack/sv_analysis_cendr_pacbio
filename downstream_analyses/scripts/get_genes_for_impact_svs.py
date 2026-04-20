from collections import defaultdict
from collections import Counter
from pathlib import Path

strain_map = {"0": "DL238", "1": "ECA36", "2" : "ECA396", "3" : "EG4725", "4" : "JU1400", "5" : "JU2526", "6" : "JU2600", "7" : "JU310", "8" : "MY2147", "9" : "MY2693", "10" : "NIC2", "11" : "NIC526", "12" : "QX1794", "13" : "XZ1516"}
IMPACT_TYPES = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

def get_strains(strains_with_sv):
	strains_with_sv = strains_with_sv.split()
	strains_with_sv = [str(i) for i, e in enumerate(strains_with_sv) if (e != "./.:NA:NA:NA:NA") and ("0/0" not in e)]
	strains_with_sv_mapped = list(map(strain_map.get, strains_with_sv, strains_with_sv))
	strains_with_sv_mapped = ', '.join(strains_with_sv_mapped)
	
	return(strains_with_sv_mapped)

impact_genes = defaultdict(list)
impact_pseudogenes = defaultdict(list)
impact_ncgenes = defaultdict(list) # Stores genes that aren't protein coding or pseudogenes
impact_genes_svnames = defaultdict(list)
impact_pseudogenes_svnames = defaultdict(list)
impact_ncgenes_svnames = defaultdict(list) # Stores genes that aren't protein coding or pseudogenes
impact_intergenic = defaultdict(list)

with open(snakemake.input[0], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t", 9)
			chromosome = line_split[0]
			start_coord = line_split[1]
			sv_name = line_split[2]
			sv_info = line_split[7]
			strains_with_sv = line_split[9]
			end_coord = None
			wormbase_ids = []
			strain_count = None
			consequences = []
			
			for info_field in sv_info.split(";"):
				if "END=" in info_field and "AVG_END=" not in info_field:
					end_coord = info_field.replace("END=", "")
				if "WBGene" in info_field:
					for x in info_field.split("|"):
						if "WBGene" in x:
							wormbase_ids.append(x)
				if "CSQ=" in info_field:
					csq_info = info_field.split("|")[1:]
					csq_info = info_field.rstrip("|").split("|")
					csq_info = info_field[:-1].split("|")
					csq_info_length = len(csq_info)

					if csq_info_length % 22 != 0:
						if (csq_info_length+2) % 22 != 0: # Entries for pseudogene, tRNA, and ncRNA seem to cause a lot of these problems
							print(line)
							print(csq_info)
							print(csq_info_length)
							input("Length off by more than 2")

					for i in range(0, csq_info_length, 22):
						end = i+22
						if end > len(csq_info):
							end = i+20
						csq_entry = csq_info[i:end]
						consequence = csq_entry[1]
						impact = csq_entry[2]
						wbgene = csq_entry[4]
						biotype = csq_entry[7]
												
						if wbgene == "":
							if consequence != "intergenic_variant":
								print("consequence: " + consequence)
								print("impact: " + impact)
								input(csq_entry)
							else:
								impact_intergenic[impact.strip()].append(sv_name)
						else:
							wbgene_svname = wbgene + "," + sv_name
							if biotype == "pseudogene":
								impact_pseudogenes[impact.strip()].append(wbgene)
								impact_pseudogenes_svnames[impact.strip()].append(wbgene_svname)
							if biotype == "protein_coding":
								impact_genes[impact.strip()].append(wbgene)
								impact_genes_svnames[impact.strip()].append(wbgene_svname)
							else:
								impact_ncgenes[impact.strip()].append(wbgene)
								impact_ncgenes_svnames[impact.strip()].append(wbgene_svname)
							
				if "SUPP=" in info_field:
					strain_count = int(info_field.split("=")[1])

			if end_coord is None:
				print("Didn't find end coordinate")

			if wormbase_ids is None:
				print("Didn't find WormBase ID")

			strains_with_sv = get_strains(strains_with_sv)
			if len(strains_with_sv.split(",")) != strain_count and strain_count > 1:
				if "0/0" not in line_split[9]:
					print("strain_count doesn't match number of strains with SV")
					print("strain_count: " + str(strain_count))
					print(strains_with_sv)
	

gene_impact_summary = ["Summary of the impact on genes/intergenic regions affected by SVs\n"]
gene_impact_summary.append("GENES\n-----\n")
for impact in impact_genes.keys():
	gene_impact_summary.append("Predicted impact: " + impact)
	gene_list = impact_genes[impact]
	gene_impact_summary.append("Total number of genes affected by this impact type: " +  str(len(gene_list)))
	gene_set = set(impact_genes[impact])
	gene_impact_summary.append("Total number of unique genes (no intergenic or pseudogenes) affected by this impact type: : " +  str(len(gene_set)) + "\n")

gene_impact_summary.append("NONCODING GENES\n-----------\n")
for impact in impact_ncgenes.keys():
	gene_impact_summary.append("Predicted impact: " + impact)
	gene_list = impact_ncgenes[impact]
	gene_impact_summary.append("Total number of noncoding genes affected by this impact type: " +  str(len(gene_list)))
	gene_set = set(impact_ncgenes[impact])
	gene_impact_summary.append("Total number of unique non coding genes affected by this impact type: : " +  str(len(gene_set)) + "\n")

gene_impact_summary.append("PSEUDOGENES\n-----------\n")
for impact in impact_pseudogenes.keys():
	gene_impact_summary.append("Predicted impact: " + impact)
	gene_list = impact_pseudogenes[impact]
	gene_impact_summary.append("Total number of pseudogenes affected by this impact type: " +  str(len(gene_list)))
	gene_set = set(impact_pseudogenes[impact])
	gene_impact_summary.append("Total number of unique pseudogenes affected by this impact type: : " +  str(len(gene_set)) + "\n")

with open(snakemake.output[0], 'w') as f:
	f.writelines(f'{s}\n' for s in gene_impact_summary)

impact_types_to_use = snakemake.params[0]
combined_impacts = '\t'.join(snakemake.params[0])

genefiles = []
ncgenefiles = []
pseudogenefiles =[]
for outfile in snakemake.output[1:]:
	if "pseudo" in outfile:
		pseudogenefiles.append(outfile)
	elif "ncgene" in outfile:
		ncgenefiles.append(outfile)	
	else:
		genefiles.append(outfile)

for impact in impact_genes.keys():
	for outfile in genefiles:
		if impact in outfile:
			if "txt" in outfile:
				gene_set = set(impact_genes[impact])
				with open(outfile, 'w') as f:
					f.writelines(f'{s}\n' for s in gene_set)
			if "csv" in outfile:
				gene_set = set(impact_genes_svnames[impact])
				with open(outfile, 'w') as f:
					f.writelines(f'{s}\n' for s in gene_set)


for impact in impact_ncgenes.keys():
	for outfile in ncgenefiles:
		if impact in outfile:
			if "txt" in outfile:
				gene_set = set(impact_ncgenes[impact])
				with open(outfile, 'w') as f:
					f.writelines(f'{s}\n' for s in gene_set)
			if "csv" in outfile:
				gene_set = set(impact_ncgenes_svnames[impact])
				with open(outfile, 'w') as f:
					f.writelines(f'{s}\n' for s in gene_set)


for impact in impact_pseudogenes.keys():
	for outfile in pseudogenefiles:
		if impact in outfile:
			if "txt" in outfile:		
				gene_set = set(impact_pseudogenes[impact])
				with open(outfile, 'w') as f:
					f.writelines(f'{s}\n' for s in gene_set)
			if "csv" in outfile:		
				gene_set = set(impact_pseudogenes_svnames[impact])
				with open(outfile, 'w') as f:
					f.writelines(f'{s}\n' for s in gene_set)

for impact in impact_types_to_use:
	if impact not in impact_genes.keys():
		print("No genes SVs for impact type: " + impact)
		for outfile in genefiles:
			if impact in outfile:
				Path(outfile).touch()
	if impact not in impact_ncgenes.keys():
		print("No noncoding genes SVs for impact type: " + impact)
		for outfile in ncgenefiles:
			if impact in outfile:
				Path(outfile).touch()	
	if impact not in impact_pseudogenes.keys():
		print("No pseudogenes SVs for impact type: " + impact)
		for outfile in pseudogenefiles:
			if impact in outfile:
				Path(outfile).touch()
