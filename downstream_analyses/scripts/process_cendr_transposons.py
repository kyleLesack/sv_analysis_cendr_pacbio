import pybedtools
from pathlib import Path

bedlines = []
with open(snakemake.input[0], 'r') as f:
	for line in f:
		line_split = line.strip().split()
		chromosome = line_split[0]
		start_coord = line_split[1]
		stop_coord = line_split[2]
		strain = line_split[3]
		transposon = line_split[4]
		
		if strain == snakemake.params.strain:
			bedline = chromosome + "\t" + start_coord + "\t" + stop_coord + "\t" + transposon
			bedlines.append(bedline)

outdir = str(Path(snakemake.output[0]).parent) 

try:
	Path(outdir).mkdir(parents=True, exist_ok=False)
except FileExistsError:
	print("Directory exists: " + outdir)
else:
	print("Directory created: " + outdir)
	
bedfile = pybedtools.BedTool(bedlines)
bedfile.sort()
outfile = snakemake.output[0]
bedfile.saveas(outfile)	
