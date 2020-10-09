###################
# count the mirtron reads from the output of contamination for small RNA-seq data (initially)
###################
DIR_bams="$1"
OUT="contamination_countTable.txt"
ml load samtools/1.4-foss-2018b

## sequence for mirtron-35 and mirtron-51
mirtron35='TCACCGGGTGTAAACTTACAG'
mirtron51='TACCCGTAATCTTCATAATTACAG'

header='sample'

for i in {20..30}; do
    header="$header mirtron35_read_${i} mirtron35_umi_${i}"
done

for i in {20..30}; do
    header="$header mirtron51_read_${i} mirtron51_umi_${i}"
done

echo $header|tr ' ' '\t' > $OUT;
echo $header

# loop over bam files
for bam in `ls ${DIR_bams}/*.bam`
do  
    bname=$(basename $bam);
    keep="$bname"
    
    # select reads mapped to mirtron-35 and of length 20-30bp, and count associated read and umi 
    mirtron35_r20=`samtools view $bam | awk '$10 == "TCACCGGGTGTAAACTTACA"' | wc -l`
    mirtron35_u20=`samtools view $bam | awk '$10 == "TCACCGGGTGTAAACTTACA"' | cut -f1|tr '_' '\t'|cut -f2|sort -u|wc -l`
    keep="$keep $mirtron35_r20 $mirtron35_u20"
    
    mirtron35_r21=`samtools view $bam | awk '$10 == "TCACCGGGTGTAAACTTACAG"' | wc -l`
    mirtron35_u21=`samtools view $bam | awk '$10 == "TCACCGGGTGTAAACTTACAG"' | cut -f1|tr '_' '\t'|cut -f2|sort -u|wc -l`
    keep="$keep $mirtron35_r21 $mirtron35_u21"
    
    for ll in {22..30};do
	#echo $ll;
	mirtron35_rr=`samtools view $bam | awk -v a=$ll '{ if(($10 ~ /^TCACCGGGTGTAAACTTACAG/) && (length($10) == a)) {print}}' | wc -l`
	mirtron35_uu=`samtools view $bam | awk -v a=$ll '{ if(($10 ~ /^TCACCGGGTGTAAACTTACAG/) && (length($10) == a)) {print}}' | cut -f1|tr '_' '\t'|cut -f2|sort -u|wc -l`

	keep="$keep $mirtron35_rr $mirtron35_uu"
    done
    
    
    # select reads mapped to mirtron-51  and of length 20-30bp, and count associated read and umi
    mirtron51_r20=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATT"' | wc -l`
    mirtron51_u20=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATT"' | cut -f1|tr '_' '\t'|cut -f2|sort -u|wc -l`
    keep="$keep $mirtron51_r20 $mirtron51_u20"
    
    mirtron51_r21=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATTA"' | wc -l`
    mirtron51_u21=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATTA"' | cut -f1|tr '_' '\t'|cut -f2|sort -u|wc -l`
    keep="$keep $mirtron51_r21 $mirtron51_u21"

    mirtron51_r22=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATTAC"' | wc -l`
    mirtron51_u22=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATTAC"' | cut -f1|tr '_' '\t'|cut -f2|sort -u|wc -l`
    keep="$keep $mirtron51_r22 $mirtron51_u22"
    
    mirtron51_r23=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATTACA"' | wc -l`
    mirtron51_u23=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATTACA"' | cut -f1|tr '_' '\t'|cut -f2|sort -u|wc -l`
    keep="$keep $mirtron51_r23 $mirtron51_u23"
    
    mirtron51_r24=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATTACAG"' | wc -l`
    mirtron51_u24=`samtools view $bam | awk '$10 == "TACCCGTAATCTTCATAATTACAG"' | cut -f1|tr '_' '\t'|cut -f2|sort -u|wc -l`
    keep="$keep $mirtron51_r24 $mirtron51_u24"
    
    for ll in {25..30};do
	#echo $ll;
	mirtron51_rr=`samtools view $bam | awk -v a=$ll '{ if(($10 ~ /^TACCCGTAATCTTCATAATTACAG/) && (length($10) == a)) {print}}' | wc -l`
	mirtron51_uu=`samtools view $bam | awk -v a=$ll '{ if(($10 ~ /^TACCCGTAATCTTCATAATTACAG/) && (length($10) == a)) {print}}' | cut -f1|tr '_' '\t'|cut -f2|sort -u|wc -l`

	keep="$keep $mirtron51_rr $mirtron51_uu"
	
    done
    
    echo $keep
    
    echo $keep|tr ' ' '\t' >> $OUT;
    
    #break
    
done
