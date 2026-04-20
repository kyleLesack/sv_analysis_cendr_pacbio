# GPCRs

grep -i "Serpentine Receptor" ../Caenorhabditis_elegans.WBcel235.61.gff3 > tmp
grep -iF "7TM_GPCR" ../Caenorhabditis_elegans.WBcel235.61.gff3 >> tmp
grep -iF "G_PROTEIN_RECEP" ../Caenorhabditis_elegans.WBcel235.61.gff3 >> tmp
sort -u tmp | grep -i "ID=gene" > gpcrs/genes_and_candidates.tsv

grep -i "Serpentine Receptor" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 > tmp
grep -iF "7TM_GPCR" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
grep -iF "G_PROTEIN_RECEP" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
sort -u tmp > gpcrs/genes_and_candidates_wb_names.txt

# F-box

grep -i "F-box" ../Caenorhabditis_elegans.WBcel235.61.gff3 > tmp
grep -iF "Name=fbxa-" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iFv "F-box" >> tmp
grep -iF "Name=fbxb-" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iFv "F-box" >> tmp
sort -u tmp | grep -i "ID=gene" > fbox/genes_and_candidates.tsv

grep -i "F-box" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 > tmp
grep -iF "Name=fbxa-" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
grep -iF "Name=fbxb-" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
sort -u tmp > fbox/genes_and_candidates_wb_names.txt

# CLEC

grep -i "C-type lectin" ../Caenorhabditis_elegans.WBcel235.61.gff3 > tmp
grep -iF "Name=clec-" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iFv "C-type lectin" >> tmp
sort -u tmp | grep -i "ID=gene" > clec/genes_and_candidates.tsv

grep -i "C-type lectin" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 > tmp
grep -iF "Name=clec-" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
sort -u tmp > clec/genes_and_candidates_wb_names.txt

# CYP

grep -i "CYtochrome P450 family" ../Caenorhabditis_elegans.WBcel235.61.gff3 > tmp
grep -i "Putative cytochrome P450" ../Caenorhabditis_elegans.WBcel235.61.gff3 >> tmp
grep -i "Cytochrome P450 daf-9" ../Caenorhabditis_elegans.WBcel235.61.gff3 >> tmp
sort -u tmp | grep -i "ID=gene" > cyp/genes_and_candidates.tsv

grep -i "CYtochrome P450 family" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 > tmp
grep -i "Putative cytochrome P450" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
grep -i "Cytochrome P450 daf-9" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
sort -u tmp > cyp/genes_and_candidates_wb_names.txt

# UGT

grep -i "ugt-" ../Caenorhabditis_elegans.WBcel235.61.gff3 > tmp
grep -i "Glucuronosyltransferase" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "3-beta-glucuronosyltransferase" | grep -iv "ugt-" >> tmp
sort -u tmp > ugt/genes_and_candidates.tsv

grep -i "ugt-" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 > tmp
grep -i "Glucuronosyltransferase" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "3-beta-glucuronosyltransferase" | grep -iv "ugt-" | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
sort -u tmp | grep -i "ID=gene" > ugt/genes_and_candidates_wb_names.txt

# GST

grep -i "gst-" ../Caenorhabditis_elegans.WBcel235.61.gff3 > tmp
grep -i "Glutathione S-Transferase" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "gst-" >> tmp
grep -i "GST" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "gst-" >> tmp
sort -u tmp | grep -i "ID=gene" > gst/genes_and_candidates.tsv

grep -i "gst-" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 > tmp
grep -i "Glutathione S-Transferase" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "gst-" | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
grep -i "GST" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "gst-" | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
sort -u tmp > gst/genes_and_candidates_wb_names.txt

# NHR

grep -i "nhr-" ../Caenorhabditis_elegans.WBcel235.61.gff3 > tmp
grep -i "Nuclear Hormone Receptor family" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "nhr-" >> tmp
grep -i "NR LBD domain-containing protein" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "nhr-" >> tmp
sort -u tmp | grep -i "ID=gene" > nhr/genes_and_candidates.tsv

grep -i "nhr-" ../Caenorhabditis_elegans.WBcel235.61.gff3 | cut -d ":" -f 2 | cut -d ";" -f 1 > tmp
grep -i "Nuclear Hormone Receptor family" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "nhr-" | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
grep -i "NR LBD domain-containing protein" ../Caenorhabditis_elegans.WBcel235.61.gff3 | grep -iv "nhr-" | cut -d ":" -f 2 | cut -d ";" -f 1 >> tmp
sort -u tmp > nhr/genes_and_candidates_wb_names.txt

rm tmp
