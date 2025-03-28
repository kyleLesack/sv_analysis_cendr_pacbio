# Create Upset plots

python scripts/3_summarize_results/create_upset.py 3_jasmine/ngmlr/sniffles/7_ins_to_dup/sv_types/DEL_SORTED.vcf 29 12 --plot
python scripts/3_summarize_results/create_upset.py 3_jasmine/ngmlr/sniffles/7_ins_to_dup/sv_types/DUP_SORTED.vcf 7 10 --plot
python scripts/3_summarize_results/create_upset.py 3_jasmine/ngmlr/sniffles/7_ins_to_dup/sv_types/INV_SORTED.vcf 5 12 --plot

# Summarize data

mkdir -p 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/
python scripts/3_summarize_results/create_upset.py 3_jasmine/ngmlr/sniffles/7_ins_to_dup/sv_types/DEL_SORTED.vcf 29 12 | grep "found" | sed 's/found/deletions/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/deletion_summary.txt
python scripts/3_summarize_results/create_upset.py 3_jasmine/ngmlr/sniffles/7_ins_to_dup/sv_types/DUP_SORTED.vcf 7 10 | grep "found" | sed 's/found/duplications/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/duplication_summary.txt
python scripts/3_summarize_results/create_upset.py 3_jasmine/ngmlr/sniffles/7_ins_to_dup/sv_types/INV_SORTED.vcf 5 12 | grep "found" | sed 's/found/inversions/g' > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/inversion_summary.txt

rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/deletion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " | sort > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/deletion_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/duplication_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " | sort > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/duplication_summary.csv
rev 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/inversion_summary.txt | cut -d " " -f 1,2 | rev | tr ":" "," | tr -d " " | sort > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/inversion_summary.csv


# Create Tables

sed -i '1i Strain, DEL' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/deletion_summary.csv
sed -i '1i Strain, DUP' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/duplication_summary.csv
sed -i '1i Strain, INV' 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/inversion_summary.csv

csvjoin --columns=1 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/deletion_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/duplication_summary.csv 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/inversion_summary.csv > 8_summary_statistics/results/ngmlr/sniffles/7_ins_to_dup/singletons/sv_summary.csv 


