# Summarize Total SV Calls

mkdir -p 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/

python3 scripts/3_summarize_results/summarize_sv_calls.py 3_jasmine/ngmlr/sniffles/7_ins_to_dup/sv_types/DEL_SORTED.vcf | sort | sed 's/variants/deletions/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/deletion_summary.txt
python3 scripts/3_summarize_results/summarize_sv_calls.py 3_jasmine/ngmlr/sniffles/7_ins_to_dup/sv_types/DUP_SORTED.vcf | sort | sed 's/variants/duplications/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/duplication_summary.txt
python3 scripts/3_summarize_results/summarize_sv_calls.py 3_jasmine/ngmlr/sniffles/7_ins_to_dup/sv_types/INV_SORTED.vcf | sort | sed 's/variants/inversions/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/inversion_summary.txt

mkdir -p 8_summary_statistics/results/ngmlr/sniffles/8_final/

python3 scripts/3_summarize_results/summarize_sv_calls.py 3_jasmine/ngmlr/sniffles/8_final/sv_types/DEL_SORTED.vcf | sort | sed 's/variants/deletions/g' > 8_summary_statistics/results/ngmlr/sniffles/8_final/deletion_summary.txt
python3 scripts/3_summarize_results/summarize_sv_calls.py 3_jasmine/ngmlr/sniffles/8_final/sv_types/DUP_SORTED.vcf | sort | sed 's/variants/duplications/g' > 8_summary_statistics/results/ngmlr/sniffles/8_final/duplication_summary.txt
python3 scripts/3_summarize_results/summarize_sv_calls.py 3_jasmine/ngmlr/sniffles/8_final/sv_types/INV_SORTED.vcf | sort | sed 's/variants/inversions/g' > 8_summary_statistics/results/ngmlr/sniffles/8_final/inversion_summary.txt

rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/deletion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/deletion_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/duplication_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/duplication_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/inversion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/inversion_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/8_final/deletion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/8_final/deletion_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/8_final/duplication_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/8_final/duplication_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/8_final/inversion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/8_final/inversion_summary.csv

# Summarize VEP High Impact SV Calls

mkdir -p 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/

python3 scripts/3_summarize_results/summarize_sv_calls.py 4_vep/output/ngmlr/7_ins_to_dup/annotated_high_impact/DEL-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf | sort | sed 's/variants/deletions/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/deletion_summary.txt
python3 scripts/3_summarize_results/summarize_sv_calls.py 4_vep/output/ngmlr/7_ins_to_dup/annotated_high_impact/DUP-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf | sort | sed 's/variants/duplications/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/duplication_summary.txt
python3 scripts/3_summarize_results/summarize_sv_calls.py 4_vep/output/ngmlr/7_ins_to_dup/annotated_high_impact/INV-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf | sort | sed 's/variants/inversions/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/inversion_summary.txt

mkdir -p 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/

python3 scripts/3_summarize_results/summarize_sv_calls.py 4_vep/output/ngmlr/8_final/annotated_high_impact/DEL-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf | sort | sed 's/variants/deletions/g' > 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/deletion_summary.txt
python3 scripts/3_summarize_results/summarize_sv_calls.py 4_vep/output/ngmlr/8_final/annotated_high_impact/DUP-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf | sort | sed 's/variants/duplications/g' > 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/duplication_summary.txt
python3 scripts/3_summarize_results/summarize_sv_calls.py 4_vep/output/ngmlr/8_final/annotated_high_impact/INV-vep_per_gene-pick-no_intergenic-annotated_high_impact.vcf | sort | sed 's/variants/inversions/g' > 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/inversion_summary.txt

rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/deletion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/deletion_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/duplication_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/duplication_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/inversion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/inversion_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/deletion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/deletion_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/duplication_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/duplication_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/inversion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " > 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/inversion_summary.csv

# Create Tables

sed -i '1i Strain, DEL' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/deletion_summary.csv
sed -i '1i Strain, DUP' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/duplication_summary.csv
sed -i '1i Strain, INV' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/inversion_summary.csv
sed -i '1i Strain, DEL' 8_summary_statistics/results/ngmlr/sniffles/8_final/deletion_summary.csv
sed -i '1i Strain, DUP' 8_summary_statistics/results/ngmlr/sniffles/8_final/duplication_summary.csv
sed -i '1i Strain, INV' 8_summary_statistics/results/ngmlr/sniffles/8_final/inversion_summary.csv
sed -i '1i Strain, DEL-VEP' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/deletion_summary.csv
sed -i '1i Strain, DUP-VEP' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/duplication_summary.csv
sed -i '1i Strain, INV-VEP' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/inversion_summary.csv
sed -i '1i Strain, DEL-VEP' 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/deletion_summary.csv
sed -i '1i Strain, DUP-VEP' 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/duplication_summary.csv
sed -i '1i Strain, INV-VEP' 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/inversion_summary.csv

csvjoin --columns=1 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/duplication_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/inversion_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/sv_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/sv_summary.csv

csvjoin --left -I --columns=1 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/duplication_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/inversion_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/vep_summary.csv

csvjoin --left -I --columns=1 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/duplication_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/duplication_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/inversion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/inversion_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/vep_high_impact/sv_vep_summary.csv

csvjoin --columns=1 8_summary_statistics/results/ngmlr/sniffles/8_final/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/8_final/duplication_summary.csv 8_summary_statistics/results/ngmlr/sniffles/8_final/inversion_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/8_final/sv_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/8_final/sv_summary.csv

csvjoin --left -I --columns=1 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/duplication_summary.csv 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/inversion_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/vep_summary.csv

csvjoin --left -I --columns=1 8_summary_statistics/results/ngmlr/sniffles/8_final/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/8_final/duplication_summary.csv 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/duplication_summary.csv 8_summary_statistics/results/ngmlr/sniffles/8_final/inversion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/inversion_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/8_final/vep_high_impact/sv_vep_summary.csv
