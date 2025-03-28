# Create Upset plots

python scripts/3_summarize_results/create_upset.py 4_vep/output/ngmlr/7_ins_to_dup/annotated_high_impact/DEL-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf 4 9 --plot
python scripts/3_summarize_results/create_upset.py 4_vep/output/ngmlr/7_ins_to_dup/annotated_high_impact/DUP-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf 2 3 --plot

# Summarize data

mkdir -p 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/
python scripts/3_summarize_results/create_upset.py 4_vep/output/ngmlr/7_ins_to_dup/annotated_high_impact/DEL-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf 4 9 | grep "found" | sed 's/found/vep deletions/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/deletion_summary.txt
python scripts/3_summarize_results/create_upset.py 4_vep/output/ngmlr/7_ins_to_dup/annotated_high_impact/DUP-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf 2 3 | grep "found" | sed 's/found/vep duplications/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/duplication_summary.txt

rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/deletion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " | sort > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/deletion_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/duplication_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " | sort > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/duplication_summary.csv

# Create Tables

sed -i '1i Strain, DEL' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/deletion_summary.csv
sed -i '1i Strain, DUP' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/duplication_summary.csv

csvjoin --columns=1 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/duplication_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/vep_high_impact/sv_summary.csv 



