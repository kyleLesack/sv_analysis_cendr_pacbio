from pathlib import Path
from collections import Counter

strain_map = {"0": "DL238", "1": "ECA36", "2" : "ECA396", "3" : "EG4725", "4" : "JU1400", "5" : "JU2526", "6" : "JU2600", "7" : "JU310", "8" : "MY2147", "9" : "MY2693", "10" : "NIC2", "11" : "NIC526", "12" : "QX1794", "13" : "XZ1516"}

def get_strains(strains_with_sv):
	strains_with_sv = strains_with_sv.split()
	strains_with_sv = [str(i) for i, e in enumerate(strains_with_sv) if (e != "./.:NA:NA:NA:NA") and ("0/0" not in e)]
	strains_with_sv_mapped = list(map(strain_map.get, strains_with_sv, strains_with_sv))
	strains_with_sv_mapped = ', '.join(strains_with_sv_mapped)
	
	return(strains_with_sv_mapped)

strain_count_list = []
with open(snakemake.input[0], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t", 9)
			start_coord = line_split[1]
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
			if len(strains_with_sv.split(",")) != strain_count:
				if "0/0" not in line_split[9]:
					print("strain_count doesn't match number of strains with SV")
					print("strain_count: " + str(strain_count))
					print(strains_with_sv)
				elif "1/0" not in line_split[9] and "0/1" not in line_split[9] and "1/1" not in line_split[9]:
					print("strain_count doesn't match number of strains with SV")
					print("strain_count: " + str(strain_count))
					print(strains_with_sv)
					input(line)
				else:
					corrected_strain_count = len(strains_with_sv.split(","))
					strain_count_list.append(corrected_strain_count)
			else:
				strain_count_list.append(strain_count)
			
			if strain_count != strain_count_ext:
				print(sv_info)
				input("SUPP != SUPP_EXT")



strain_count_list_counter = Counter(strain_count_list)
for x in range(1,15):
	if x not in strain_count_list_counter.keys():
		strain_count_list_counter[x] = 0

strain_count_list_counter = dict(sorted(strain_count_list_counter.items()))

strain_count_csv = ["Strains,Count"]
for strain_count in strain_count_list_counter:
	csv_line = str(strain_count) + "," + str(strain_count_list_counter[strain_count])
	strain_count_csv.append(csv_line)

with open(snakemake.output[0], 'w') as f:
	f.writelines(f'{s}\n' for s in strain_count_csv)