###################
# count the mirtron reads from the output of contamination for small RNA-seq data (initially)
###################
DIR_bams="$1"
OUT="mirtron_tails_countTable.txt"
ml load samtools/1.4-foss-2018b

## sequence for mirtron-35 and mirtron-51
mirtron35='TCACCGGGTGTAAACTTACAG'
mirtron51='TACCCGTAATCTTCATAATTACAG'

header='sample A C G T N'

echo $header|tr ' ' '\t' > $OUT;
echo $header

# loop over bam files
for bam in `ls ${DIR_bams}/*.bam`
do  
    bname=$(basename $bam);
        
    # select reads mapped to mirtron-35 and of length 20-30bp, and count associated read and umi 
    tail_A=`samtools view $bam |awk -v a=22 '{ if(($10 ~ /^TCACCGGGTGTAAACTTACAG/) && (length($10) == a)) {print}}'|cut -f10|grep TCACCGGGTGTAAACTTACAGA|wc -l`
    tail_C=`samtools view $bam |awk -v a=22 '{ if(($10 ~ /^TCACCGGGTGTAAACTTACAG/) && (length($10) == a)) {print}}'|cut -f10|grep TCACCGGGTGTAAACTTACAGC|wc -l`
    tail_G=`samtools view $bam |awk -v a=22 '{ if(($10 ~ /^TCACCGGGTGTAAACTTACAG/) && (length($10) == a)) {print}}'|cut -f10|grep TCACCGGGTGTAAACTTACAGG|wc -l`
    tail_T=`samtools view $bam|awk -v a=22 '{ if(($10 ~ /^TCACCGGGTGTAAACTTACAG/) && (length($10) == a)) {print}}'|cut -f10|grep TCACCGGGTGTAAACTTACAGT|wc -l`
    tail_N=`samtools view $bam|awk -v a=22 '{ if(($10 ~ /^TCACCGGGTGTAAACTTACAG/) && (length($10) == a)) {print}}'|cut -f10|grep TCACCGGGTGTAAACTTACAGN|wc -l`
    
    keep="$bname $tail_A $tail_C $tail_G $tail_T $tail_N"
    echo $keep
    echo $keep|tr ' ' '\t' >> $OUT;
    
    #break
    
done
