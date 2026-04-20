from collections import defaultdict
from pathlib import Path

strain_map = {"0": "DL238", "1": "ECA36", "2" : "ECA396", "3" : "EG4725", "4" : "JU1400", "5" : "JU2526", "6" : "JU2600", "7" : "JU310", "8" : "MY2147", "9" : "MY2693", "10" : "NIC2", "11" : "NIC526", "12" : "QX1794", "13" : "XZ1516"}

IMPACT_TYPES = ["HIGH", "MODERATE", "LOW", "MODIFIER"]


# Get list of strains with a given structural variant
def get_strains(strains_with_sv):
	strains_with_sv = strains_with_sv.split()
	strains_with_sv = [str(i) for i, e in enumerate(strains_with_sv) if (e != "./.:NA:NA:NA:NA") and ("0/0" not in e)]
	strains_with_sv_mapped = list(map(strain_map.get, strains_with_sv, strains_with_sv))
	strains_with_sv_mapped = ','.join(strains_with_sv_mapped)
	
	return(strains_with_sv_mapped)


vcf_lines = []
svs_homozygous_reference_only = 0

# Parse VCF file line by line and process variant calls
with open(snakemake.input[0], 'r') as f:
	for line in f:
		if snakemake.params[0] in line:
			print("Excluding line due to invalid genotypes:")
			print(line)
		elif "#" in line:
			vcf_lines.append(line)
		else:
			line_split = line.split("\t", 9)
			strains_with_sv = line_split[9]
			strains_with_sv = get_strains(strains_with_sv)
			if len(strains_with_sv) > 0:
				vcf_lines.append(line)
			else:
				svs_homozygous_reference_only += 1
				

outdir = str(Path(snakemake.output[0]).parent)

try:
	Path(outdir).mkdir(parents=True, exist_ok=False)
except FileExistsError:
	print("Directory exists: " + outdir)
else:
	print("Directory created: " + outdir)

outfile = Path(snakemake.output[0])

print("Number of SVs with only homozygous reference calls: " + str(svs_homozygous_reference_only))

with open(outfile, 'w') as f:
	f.writelines(s for s in vcf_lines)
