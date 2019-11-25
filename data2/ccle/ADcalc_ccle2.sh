#!/bin/bash
# ADcalc_ccle2.sh
# Allelic Depth Calculator
# for a position
#
# s3://crownproject/ccle/scripts/ADcalc_ccle2.sh

# Type [rna | wgs]
TYPE='ccle.rna'

# Controls -----------------
DEPTH='100000' #Max per file DP

# Regions in hgr1.fa reference genome
REGIONS=('chr13:1003660-1005529' 'chr13:1005529-1005629' \
        'chr13:10219-10340' 'chr13:1006622-1006779' 'chr13:1007948-1013018')

# Corresponding region/gene names
GENES=('18S' '18SE' '5S' '5.8S' '28S')

# 18S  1870
# 18SE 101
# 5S   122
# 5.8S 158
# 28S  5071

# Terminate instances upon completion (for debugging)
TERMINATE='FALSE'

# S3 Output directory
S3DIR='s3://crownproject/ccle/gvcf/'

# Script ------------------
BAMLIST='bam.list.tmp'

cd ~/ccle/
mkdir -p GVCF #Output Folder

cd hgr1

#for TYPE in $(echo "hgr1")
#do
    echo Analyzing $TYPE...
    #cd $TYPE

    ls *.wgs.hgr1.bam > bam.list.tmp
    ls *.wgs.hgr1.bam > ../GVCF/$TYPE.bamlist
          
    for index in ${!GENES[*]}
    do
      printf "Started processing %s\n" ${GENES[$index]}
      OUTPUT="../GVCF/$TYPE.${GENES[$index]}.gvcf"

      # Iterate through every bam file in directory
      # look-up position and return VCF
      bcftools mpileup -f ~/resources/hgr1/hgr1.fa \
        --max-depth $DEPTH -A --min-BQ 30 \
        -a FORMAT/DP,AD \
        -r ${REGIONS[$index]} \
        --ignore-RG \
        -b $BAMLIST | \
        bcftools annotate -x INFO,FORMAT/PL - | \
        bcftools view -O v - \
        >> $OUTPUT

      RESULTS+=("$OUTPUT")
      printf "Done with %s \n" ${GENES[$index]}
      printf "%s\n" ${REGIONS[$index]}

    done

    rm bam.list.tmp

#    cd .. # move to tcga folder to reset
#done

# Copy GVCF output to AWS S3
cd ../GVCF
aws s3 cp --recursive ./ $S3DIR

# Shutdown and Terminate instance
EC2ID=$(ec2metadata --instance-id)
sleep 20s # to catch errors

if [ "$TERMINATE" = TRUE ]
then
  echo "Run Complete -- Shutting down instance."
  aws ec2 terminate-instances --instance-ids $EC2ID
else
  echo "Run Complete -- Instance is online."
fi
