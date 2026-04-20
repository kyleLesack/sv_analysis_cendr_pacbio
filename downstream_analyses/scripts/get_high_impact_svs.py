high_impact_svs_set = set()
with open(snakemake.input[0], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t", 9)
			chromosome = line_split[0]
			start_coord = line_split[1]
			sv_name = line_split[2]
			sv_info = line_split[7]
			
			for info_field in sv_info.split(";"):
				if "HIGH" in info_field:		
					high_impact_svs_set.add(sv_name)

with open(snakemake.output[0], 'w') as f:
	f.writelines(f'{s}\n' for s in high_impact_svs_set)

