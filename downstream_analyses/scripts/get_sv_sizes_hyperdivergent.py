from pathlib import Path

hyperdivergent_svs = []
with open(snakemake.input["bedfile"], 'r') as f:
	for line in f:
		sv_name = line.split("\t")[3]
		hyperdivergent_svs.append(sv_name)
	
hyperdivergent_svs = set(hyperdivergent_svs)
vcf_size_hyperdivergent = []
vcf_size_not_hyperdivergent = []
with open(snakemake.input["vcffile"], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t", 9)
			sv_name = line_split[2]
			sv_info = line_split[7]
			sv_length = None
			
			for info_field in sv_info.split(";"):
				if "SVLEN=" in info_field:
					sv_length = abs(int(info_field.split("=")[1]))
			if sv_length is not None:
				new_line = sv_name + "\t" + str(sv_length)
				if sv_name in hyperdivergent_svs:
					vcf_size_hyperdivergent.append(new_line)
				else:
					vcf_size_not_hyperdivergent.append(new_line)
			else:
				print("Size not found")

outfile = snakemake.output["hyperfile"]
with open(outfile, 'w') as f:
	f.writelines(f'{s}\n' for s in vcf_size_hyperdivergent)
	
outfile = snakemake.output["nonhyperfile"]
with open(outfile, 'w') as f:
	f.writelines(f'{s}\n' for s in vcf_size_not_hyperdivergent)
	