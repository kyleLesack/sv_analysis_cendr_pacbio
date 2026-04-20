from collections import defaultdict
from collections import Counter
from pathlib import Path

strain_map = {"0": "DL238", "1": "ECA36", "2" : "ECA396", "3" : "EG4725", "4" : "JU1400", "5" : "JU2526", "6" : "JU2600", "7" : "JU310", "8" : "MY2147", "9" : "MY2693", "10" : "NIC2", "11" : "NIC526", "12" : "QX1794", "13" : "XZ1516"}

IMPACT_TYPES = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

with open(snakemake.input["referencefile"], 'r') as f:
	header = f.readline()
	vep_consequences_impacts = f.readlines()

impact_consequences_dict = defaultdict(list)
for line in vep_consequences_impacts:
	splitline = line.split("\t")
	consequence = splitline[0].strip()
	impact = splitline[4].strip()
	impact_consequences_dict[impact].append(consequence)

def get_consequences_for_impact(predicted_consequences, impact_consequences_dict):
	predicted_consequences_dict = {}
	for impact in impact_consequences_dict.keys():
		vep_consequences = impact_consequences_dict[impact]
		impact_consequences = [value for value in predicted_consequences if value in vep_consequences]
		predicted_consequences_dict[impact] = impact_consequences

	return(predicted_consequences_dict)

def get_strains(strains_with_sv):
	strains_with_sv = strains_with_sv.split()
	strains_with_sv = [str(i) for i, e in enumerate(strains_with_sv) if (e != "./.:NA:NA:NA:NA") and ("0/0" not in e)]
	strains_with_sv_mapped = list(map(strain_map.get, strains_with_sv, strains_with_sv))
	strains_with_sv_mapped = ', '.join(strains_with_sv_mapped)
	
	return(strains_with_sv_mapped)

impact_consequences = defaultdict(list)
with open(snakemake.input["vcffile"], 'r') as f:
	for line in f:
		if "WBGene" in line or "intergenic" in line:
			line_split = line.split("\t", 9)
			chromosome = line_split[0]
			start_coord = line_split[1]
			sv_name = line_split[2]
			sv_info = line_split[7]
			strains_with_sv = line_split[9]
			end_coord = None
			wormbase_ids = []
			predicted_impact = []
			strain_count = None
			
			for info_field in sv_info.split(";"):
				if "END=" in info_field and "AVG_END=" not in info_field:
					end_coord = info_field.replace("END=", "")
				if "WBGene" in info_field:
					for x in info_field.split("|"):
						if "WBGene" in x:
							wormbase_ids.append(x)
				if "CSQ=" in info_field:
					#csq_info = print(line)
					csq_info = info_field.split("|")[1:]
					#print(len(info_field.split("|")))
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
						impact = csq_entry[2].strip()
						consequence_split = consequence.split("&")

						if "&" in consequence:
							predicted_consequences_dict = get_consequences_for_impact(consequence_split, impact_consequences_dict)

							for impact in predicted_consequences_dict.keys():
								consequences = predicted_consequences_dict[impact]
								if (len(consequences)) > 0:
									for x in consequences:
										impact_consequences[impact].append(x)
						else:
							impact_consequences[impact].append(consequence)
						
					for x in info_field.split("|"):
						if x in IMPACT_TYPES:
							predicted_impact.append(x.strip())
							
			if end_coord is None:
				print("Didn't find end coordinate")

csv_lines = ["Severity,Consequence,Count"]
high_impact_consequences = impact_consequences["HIGH"]
consequence_count = Counter(impact_consequences["HIGH"]).most_common()
for line in consequence_count:
	consequence = line[0]
	count = line[1]
	csv_line = "HIGH," + consequence + "," + str(count)
	csv_lines.append(csv_line)

moderate_impact_consequences = impact_consequences["MODERATE"]
consequence_count = Counter(impact_consequences["MODERATE"]).most_common()
for line in consequence_count:
	consequence = line[0]
	count = line[1]
	csv_line = "MODERATE," + consequence + "," + str(count)
	csv_lines.append(csv_line)

low_impact_consequences = impact_consequences["LOW"]
consequence_count = Counter(impact_consequences["LOW"]).most_common()
for line in consequence_count:
	consequence = line[0]
	count = line[1]
	csv_line = "LOW," + consequence + "," + str(count)
	csv_lines.append(csv_line)

modifier_impact_consequences = impact_consequences["MODIFIER"]
consequence_count = Counter(impact_consequences["MODIFIER"]).most_common()
for line in consequence_count:
	consequence = line[0]
	count = line[1]
	csv_line = "MODIFIER," + consequence + "," + str(count)
	csv_lines.append(csv_line)

with open(snakemake.output[0], 'w') as f:
	f.writelines(f'{s}\n' for s in csv_lines)
