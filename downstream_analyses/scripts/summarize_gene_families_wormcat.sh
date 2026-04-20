genefamily=$1
vepgenes=$2
wormcatfile=$3
outdir=5_wormcat_ora/gene_families/$genefamily

wbgenefile=$outdir/wbgenes.txt
summary_file=$outdir/summary.txt
cat $vepgenes | cut -d "," -f 1 | sort -u > $wbgenefile
echo -n "Total unique $genefamily genes: " > $summary_file
cat $wbgenefile | wc -l >> $summary_file

IFS=$'\n'
for i in $(grep -iFf $wbgenefile $wormcatfile | cut -d "," -f 1 | tr -d '"'); do
   echo -n "$i: " >> $summary_file
   grep "$i" $wormcatfile | grep -ioFf $wbgenefile | wc -l >> $summary_file
   wormcat=${i// /_}   
   wormcat="${wormcat%"${wormcat##*[![:space:]]}"}" 
   wcfile=$outdir/${wormcat}_$genefamily.txt
   grep "$i" $wormcatfile | grep -ioFf $wbgenefile > $wcfile
   wcfileinfo=$outdir/${wormcat}_${genefamily}_svinfo.txt   
   grep -iFf $wcfile $vepgenes > $wcfileinfo
done