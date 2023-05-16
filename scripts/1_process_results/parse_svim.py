import pysam
import os

# Extract svim variants
def parse_svim(svim_variants):
	vcf_header_line = []
	vcf_lines = []
	excluded_due_to_type = []
	excluded_due_to_size = []
	excluded_due_to_filter = []
	excluded_due_to_qual = []

	for line in svim_variants:
		if line[0] == "#":
			if "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample" in line:
				line = line.replace("Sample", snakemake.params['sample_name'])
				vcf_header_line.append(line)
			else:
				vcf_header_line.append(line)

		else:
			line_split = line.split()
			chromosome = line_split[0]
			start_coord = line_split[1]
			#variant_name = line_split[2]
			#line = line.replace(variant_name, "svim_sv")
			variant_qual = int(line_split[5]) # svim qual score
			filter_value = line_split[6] # Describes if variant call passed. Most calls are designated as PASSED, so I probably can ignore the rest, which are probably low quality anyways.
			info_line = line_split[7] # Get info field metadata
			info_line_split = info_line.split(";")
			end_coord = None
			variant_type = None
			variant_size = None

			if filter_value != "PASS":
				excluded_due_to_filter.append(line)
			else:
				# Get the end coordinate, variant type, and size
				for x in info_line_split:
					if "END=" in x:
						end_coord = x.split("=")[1]
						variant_size = int(end_coord) - int(start_coord)
					elif "SVTYPE=" in x:
						variant_type = x.split("=")[1]
					elif "SVLEN=" in x:
						variant_size_vcf = x.split("=")[1]

				if variant_type is not None:
					if variant_type not in snakemake.params['SVIM_SV_TYPES']:
						excluded_due_to_type.append(line)
					else:
						variant_size_from_coords = int(end_coord) - int(start_coord)
						abs_variant_size= abs(int(variant_size_from_coords))
						if variant_qual >= int(snakemake.params['min_qual_svim']) and abs_variant_size >= int(snakemake.params['minsize']):
							vcf_lines.append(line)
						elif variant_qual < int(snakemake.params['min_qual_svim']):
							excluded_due_to_qual.append(line)
						else:
							excluded_due_to_size.append(line)
				else:
					print("Variant type is None")
					print(line.rstrip("\n"))

	return (vcf_header_line, vcf_lines,excluded_due_to_filter, excluded_due_to_type, excluded_due_to_size, excluded_due_to_qual)

# Write svs that were excluded to disk
def write_vcf(vcf_lines, outdir, filename):
	# Check if output directories exist.
	outdir = os.path.dirname(outdir)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	output_file = outdir + "/" + filename
	print("Writing to: " + output_file)
	with open(output_file, 'w', newline='\n') as f:
		f.writelines("%s\n" % l.rstrip() for l in vcf_lines)
	if filename == "passed.sorted.vcf":
		output_file_compressed = output_file + ".gz"
		print("Writing to " + output_file_compressed)
		pysam.tabix_compress(output_file, output_file_compressed, force=False)
		output_file_indexed = output_file_compressed + ".tbi"
		print("Writing to " + output_file_indexed)
		pysam.tabix_index(output_file_compressed, preset="vcf", force=True)

with open(snakemake.input[0]) as f:
	svim_variants = f.readlines()
	svim_variants_parsed = parse_svim(svim_variants)
	svim_header = svim_variants_parsed[0]
	svim_passed = svim_variants_parsed[1]
	passed_vcf_lines = svim_header + svim_passed
	outdir = os.path.dirname(snakemake.params['outdir'])
	passed_dir = outdir + "/passed/"
	write_vcf(passed_vcf_lines, passed_dir, "passed.sorted.vcf")
	svim_excluded_due_to_filter = svim_variants_parsed[2]
	svim_excluded_due_to_type = svim_variants_parsed[3]
	svim_excluded_due_to_size = svim_variants_parsed[4]
	svim_excluded_due_to_qual = svim_variants_parsed[5]

	excluded_due_to_filter_vcf_lines = svim_header + svim_excluded_due_to_filter
	excluded_filter_dir = outdir + "/excluded/filter/"
	write_vcf(excluded_due_to_filter_vcf_lines, excluded_filter_dir, "excluded_filter_variants.vcf")
	excluded_due_to_type_vcf_lines = svim_header + svim_excluded_due_to_type
	excluded_type_dir = outdir + "/excluded/type/"
	write_vcf(excluded_due_to_type_vcf_lines, excluded_type_dir, "excluded_type_variants.vcf")
	excluded_due_to_size_vcf_lines = svim_header + svim_excluded_due_to_size
	excluded_size_dir = outdir + "/excluded/size/"
	write_vcf(excluded_due_to_size_vcf_lines, excluded_size_dir, "excluded_size_variants.vcf")
	excluded_due_to_qual_vcf_lines = svim_header + svim_excluded_due_to_qual
	excluded_qual_dir = outdir + "/excluded/qual/"
	write_vcf(excluded_due_to_qual_vcf_lines, excluded_qual_dir, "excluded_qual_variants.vcf")
