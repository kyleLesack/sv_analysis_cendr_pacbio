import pybedtools
from collections import defaultdict
from pathlib import Path

strain_map = {"0": "DL238", "1": "ECA36", "2" : "ECA396", "3" : "EG4725", "4" : "JU1400", "5" : "JU2526", "6" : "JU2600", "7" : "JU310", "8" : "MY2147", "9" : "MY2693", "10" : "NIC2", "11" : "NIC526", "12" : "QX1794", "13" : "XZ1516"}

# Get strains with a given structural variant
def get_strains(strains_with_sv):
	strains_with_sv = strains_with_sv.split()
	strains_with_sv = [str(i) for i, e in enumerate(strains_with_sv) if (e != "./.:NA:NA:NA:NA") and ("0/0" not in e)]
	strains_with_sv_mapped = list(map(strain_map.get, strains_with_sv, strains_with_sv))
	strains_with_sv_mapped = ','.join(strains_with_sv_mapped)
	
	return(strains_with_sv_mapped)

strain_beds = defaultdict(list)
with open(snakemake.input[0], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t", 9)
			chromosome = line_split[0]
			start_coord = int(line_split[1])
			sv_name = line_split[2]
			sv_info = line_split[7]
			strains_with_sv = line_split[9]
			end_coord = None
			wormbase_ids = []
			
			# Get end coordinate and WormBase id from info field
			for info_field in sv_info.split(";"):
				if "END=" in info_field and "AVG_END=" not in info_field:
					end_coord = info_field.replace("END=", "")
				if "WBGene" in info_field:
					for x in info_field.split("|"):
						if "WBGene" in x:
							wormbase_ids.append(x)

				if "SUPP=" in info_field:
					strain_count = int(info_field.split("=")[1])

			if end_coord is None:
				print("Didn't find end coordinate")

			if wormbase_ids is None:
				print("Didn't find WormBase ID")

			strains_with_sv = get_strains(strains_with_sv)

			for strain in strains_with_sv.split(","):
				bedline = chromosome + "\t" + str(start_coord - 1) + "\t" + end_coord + "\t" + sv_name
				strain_beds[strain.strip()].append(bedline)
			
outdir = str(Path(snakemake.output[0]).parent)

try:
	Path(outdir).mkdir(parents=True, exist_ok=False)
except FileExistsError:
	print("Directory exists: " + outdir)
else:
	print("Directory created: " + outdir)

if "DEL" in snakemake.input[0]:
	svtype = "DEL"
elif "DUP" in snakemake.input[0]:
	svtype = "DUP"
elif "INV" in snakemake.input[0]:
	svtype = "INV"
else:
	svtype = "UNKNOWN"

for strain in strain_beds.keys():
	bedlines = strain_beds[strain]
	bedfile = pybedtools.BedTool(bedlines)
	bedfile.sort()
	outfile = outdir + "/" + strain + "_" + svtype + "-no_homozygous_reference.bed"
	bedfile.saveas(outfile)
