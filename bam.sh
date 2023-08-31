sed '/chrM/d;/random/d;/chrUn/d' $1_filtered.sam
samtools view -S -b $1_filtered.sam > $1.bam
samtools view -c $1.bam