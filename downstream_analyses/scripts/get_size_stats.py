from pathlib import Path
from statistics import mean, median, stdev

deletion_sizes = []
max_deletion_size = None
min_deletion_size = None
with open(snakemake.input[0], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t", 9)
			start_coord = line_split[1]
			end_coord = None
			sv_length = None
			sv_info = line_split[7]

			for info_field in sv_info.split(";"):
				if "END=" in info_field and "AVG_END=" not in info_field:
					end_coord = info_field.replace("END=", "")
				elif "SVLEN=" in info_field:
					sv_length = abs(int(info_field.split("=")[1]))

			# Calculate SV length from coordinates
			length_from_coords = int(end_coord) - int(start_coord)
			if length_from_coords != sv_length:
				print(line)
				input("SVLEN doesn't equal length calculated from coordinates")
			else:
				deletion_sizes.append(sv_length)
				
				if min_deletion_size is None:
					min_deletion_size = sv_length
				elif sv_length < min_deletion_size:
					min_deletion_size = sv_length
					
				if max_deletion_size is None:
					max_deletion_size = sv_length
				elif sv_length > max_deletion_size:
					max_deletion_size = sv_length
				
			if end_coord is None:
				print("Didn't find end coordinate")

duplication_sizes = []
max_duplication_size = None
min_duplication_size = None
with open(snakemake.input[1], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t", 9)
			start_coord = line_split[1]
			end_coord = None
			sv_length = None
			sv_info = line_split[7]

			for info_field in sv_info.split(";"):
				if "END=" in info_field and "AVG_END=" not in info_field:
					end_coord = info_field.replace("END=", "")
				elif "SVLEN=" in info_field:
					sv_length = abs(int(info_field.split("=")[1]))


			length_from_coords = int(end_coord) - int(start_coord)
			if length_from_coords != sv_length:
				print(line)
				input("SVLEN doesn't equal length calculated from coordinates")
			else:
				duplication_sizes.append(sv_length)
				
				if min_duplication_size is None:
					min_duplication_size = sv_length
				elif sv_length < min_duplication_size:
					min_duplication_size = sv_length
					
				if max_duplication_size is None:
					max_duplication_size = sv_length
				elif sv_length > max_duplication_size:
					max_duplication_size = sv_length
				
			if end_coord is None:
				print("Didn't find end coordinate")

inversion_sizes = []
max_inversion_size = None
min_inversion_size = None
with open(snakemake.input[2], 'r') as f:
	for line in f:
		if "#" not in line:
			line_split = line.split("\t", 9)
			start_coord = line_split[1]
			end_coord = None
			sv_length = None
			sv_info = line_split[7]

			for info_field in sv_info.split(";"):
				if "END=" in info_field and "AVG_END=" not in info_field:
					end_coord = info_field.replace("END=", "")
				elif "SVLEN=" in info_field:
					sv_length = abs(int(info_field.split("=")[1]))


			length_from_coords = int(end_coord) - int(start_coord)
			if length_from_coords != sv_length:
				print(line)
				input("SVLEN doesn't equal length calculated from coordinates")
			else:
				inversion_sizes.append(sv_length)
				
				if min_inversion_size is None:
					min_inversion_size = sv_length
				elif sv_length < min_inversion_size:
					min_inversion_size = sv_length
					
				if max_inversion_size is None:
					max_inversion_size = sv_length
				elif sv_length > max_inversion_size:
					max_inversion_size = sv_length
				
				
			if end_coord is None:
				print("Didn't find end coordinate")

summary_lines = ["SVType,Count,Mean_sd,Min_size,Max_size"]
summary_lines.append("DEL," + str(len(deletion_sizes)) + "," + str(round(mean(deletion_sizes))) + "±" + str(round(stdev(deletion_sizes))) + "," + str(min_deletion_size) + "," + str(max_deletion_size))
summary_lines.append("DUP," + str(len(duplication_sizes)) + "," + str(round(mean(duplication_sizes))) + "±" + str(round(stdev(duplication_sizes))) + "," + str(min_duplication_size) + "," + str(max_duplication_size))
summary_lines.append("INV," + str(len(inversion_sizes)) + "," + str(round(mean(inversion_sizes))) + "±" + str(round(stdev(inversion_sizes))) + "," + str(min_inversion_size) + "," + str(max_inversion_size))

outdir = str(Path(snakemake.output[0]).parent)

try:
	Path(outdir).mkdir(parents=True, exist_ok=False)
except FileExistsError:
	print("Directory exists: " + outdir)
else:
	print("Directory created: " + outdir)

outfile = Path(snakemake.output[0])
with open(outfile, 'w') as f:
	f.writelines(line + '\n' for line in summary_lines)