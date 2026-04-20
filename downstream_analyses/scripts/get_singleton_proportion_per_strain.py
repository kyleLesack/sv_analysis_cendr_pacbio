import argparse
from collections import defaultdict
from pathlib import Path
from collections import Counter

strain_map = {"0": "DL238", "1": "ECA36", "2" : "ECA396", "3" : "EG4725", "4" : "JU1400", "5" : "JU2526", "6" : "JU2600", "7" : "JU310", "8" : "MY2147", "9" : "MY2693", "10" : "NIC2", "11" : "NIC526", "12" : "QX1794", "13" : "XZ1516"}
STRAINS = ["DL238", "ECA36", "ECA396", "EG4725", "JU1400", "JU2526", "JU2600", "JU310", "MY2147", "MY2693", "NIC2", "NIC526", "QX1794", "XZ1516"]
HAWAIIAN_STRAINS = ["DL238", "ECA396","QX1794", "XZ1516"]
EUROPEAN_STRAINS = ["EG4725", "JU1400", "JU2526", "JU2600", "JU310", "MY2147", "MY2693", "NIC2", "NIC526"]
NEW_ZEALAND_STRAINS = ["ECA36"]

def get_strains(strains_with_sv):
	strains_with_sv = strains_with_sv.split()
	strains_with_sv = [str(i) for i, e in enumerate(strains_with_sv) if (e != "./.:NA:NA:NA:NA") and ("0/0" not in e)]
	strains_with_sv_mapped = list(map(strain_map.get, strains_with_sv, strains_with_sv))
	strains_with_sv_mapped = ','.join(strains_with_sv_mapped)
	return(strains_with_sv_mapped)

strain_svs = defaultdict(list)
svs_in_all_strains = []
with open(snakemake.input[0], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t", 9)
			start_coord = line_split[1]
			sv_name = line_split[2]
			end_coord = None
			sv_length = None
			sv_info = line_split[7]
			strains_with_sv = line_split[9]
			strain_count = None
			strain_count_ext = None

			for info_field in sv_info.split(";"):
				if "END=" in info_field and "AVG_END=" not in info_field:
					end_coord = info_field.replace("END=", "")
				elif "SVLEN=" in info_field:
					sv_length = abs(int(info_field.split("=")[1]))
				elif "SUPP=" in info_field:
					strain_count = abs(int(info_field.split("=")[1]))
				elif "SUPP_EXT=" in info_field:
					strain_count_ext = abs(int(info_field.split("=")[1]))
				
			length_from_coords = int(end_coord) - int(start_coord)
			if length_from_coords != sv_length:
				print(line)
				input("SVLEN doesn't equal length calculated from coordinates")
				
			if end_coord is None:
				print("Didn't find end coordinate")

			strains_with_sv = get_strains(strains_with_sv)
			if set(STRAINS) == set(strains_with_sv.split(",")):
				svs_in_all_strains.append(sv_name)
			
			if len(strains_with_sv.split(",")) != strain_count:
				if "0/0" not in line_split[9]:
					print("strain_count doesn't match number of strains with SV")
					print("strain_count: " + str(strain_count))
					print(strains_with_sv)
					input("Press any key to continue")
			for strain in strains_with_sv.split(","):
				strain_svs[strain.strip()].append(strain_count)
		
			if strain_count != strain_count_ext:
				print(sv_info)
				input("SUPP != SUPP_EXT")
		
csv_lines = []
for strain in STRAINS:
	singleton_count = 0
	strain_counts_for_strain =  strain_svs[strain]
		
	for count in strain_counts_for_strain:
		if count == 1:
			singleton_count += 1
	
	total_svs = len(strain_counts_for_strain)
	singleton_prop = round(singleton_count/total_svs, 2)
	if strain in HAWAIIAN_STRAINS:
		region = "Hawaii"
	elif strain in EUROPEAN_STRAINS:
		region = "Western Europe"
	elif strain in NEW_ZEALAND_STRAINS:
		region = "New Zealand"
	else:
		print("Strain region not found for " + strain) 

	csv_lines.append(region + "," + strain + "," + str(total_svs) + "," + str(singleton_count) + "," + str(singleton_prop))
		
csv_lines.sort()
with open(snakemake.output["csvfile"], 'w') as f:
	f.write("Region,Strain,Total,Singletons,Singleton_prop\n")
	f.writelines(f'{s}\n' for s in csv_lines)

with open(snakemake.output["textfile"], 'w') as f:
	f.writelines(f'{s}\n' for s in svs_in_all_strains)
	
