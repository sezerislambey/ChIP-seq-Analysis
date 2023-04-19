# Download sratoolkit:
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

# Extract sra toolkit:
tar -zxvf sratoolkit.tar.gz

# install bowtie2:
sudo apt install bowtie2

# install bedtools
sudo apt-get install bedtools

# Install MACS3 to user directory:
pip install MACS3
macs3 callpeak

# install ceas:
git clone https://github.com/taoliu/CEASp3.git
cd CEASp3
python setup.py build
pip install .

# download homer:
# http://homer.ucsd.edu/homer/ 
# Make sure that the required software specified on the site is installed on your system.
mkdir homer
mv configureHomer.pl homer
perl homer/configureHomer.pl -install
perl homer/.//configureHomer.pl -install hg19

# Download bedGraphToBigWig script:
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes

# Generate fastq files: (for example)
sratoolkit.3.0.1-ubuntu64/bin/fastq-dump -Z SRR502329 > STAT1_30m_IFNa.fastq
sratoolkit.3.0.1-ubuntu64/bin/fastq-dump -Z SRR502327 > STAT1_6h_IFNa.fastq
sratoolkit.3.0.1-ubuntu64/bin/fastq-dump -Z SRR502228 > INP_30m_IFNa.fastq
sratoolkit.3.0.1-ubuntu64/bin/fastq-dump -Z SRR502225 > INP_6h_IFNa.fastq

# Download hg19 bowtie2 index files:
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip

# Unzip hg19 bowtie2 index files:
unzip hg19.zip

# Align each fastq file to hg19 refernce:
bowtie2 -q -p 4 -k 1 --local --no-unal -x hg19 STAT1_30m_IFNa.fastq > STAT1_30m_IFNa.sam
bowtie2 -q -p 4 -k 1 --local --no-unal -x hg19 STAT1_6h_IFNa.fastq > STAT1_6h_IFNa.sam
bowtie2 -q -p 4 -k 1 --local --no-unal -x hg19 INP_30m_IFNa.fastq > INP_30m_IFNa.sam
bowtie2 -q -p 4 -k 1 --local --no-unal -x hg19 INP_6h_IFNa.fastq > INP_6h_IFNa.sam

"""
# Change permissions on file to be executable:
#chmod 755 normalize_sam.sh

# Run script on each sam file to sample 1.1x107 ChIP reads and 2x107 Input reads:
#./normalize_sam.sh STAT1_30m_IFNa.sam 11000000
#./normalize_sam.sh STAT1_6h_IFNa.sam 11000000
#./normalize_sam.sh INP_30m_IFNa.sam 19000000
#./normalize_sam.sh INP_6h_IFNa.sam 19000000
"""

# Change permissions and run shell script:
chmod 755 bam.sh
./bam.sh STAT1_30m_IFNa
./bam.sh STAT1_6h_IFNa
./bam.sh INP_30m_IFNa	# 11000000/19048444 = 0.577
./bam.sh INP_6h_IFNa	# 11000000/19670300 = 0.559

# Sub sample reads to obtain -11 million mapped reads for each input sample:
samtools view -b -s 1.577 INP_30m_IFNa.bam > INP_30m_IFNa_11E6.bam
samtools view -b -s 1.559 INP_6h_IFNa.bam > INP_6h_IFNa_11E6.bam
samtools view -c INP_6h_IFNa_11E6.bam

chmod 777 INP_*11E6.bam
chmod 777 STAT1_30m_IFNa.bam
chmod 777 STAT1_6h_IFNa.bam

# MACS3 peak calling:
# macs3 callpeak -t norm_STAT1_30m_IFNa.sam -c norm_INP_30m_IFNa.sam -n STAT1_30m_IFNa -g hs --bdg -q 0.05 -f SAM
# macs3 callpeak -t norm_STAT1_6h_IFNa.sam -c norm_INP_6h_IFNa.sam -n STAT1_6h_IFNa -g hs --bdg -q 0.05 -f SAM
macs3 callpeak -t STAT1_6h_IFNa.bam -c INP_6h_IFNa_11E6.bam -n STAT1_6h_IFNa -g hs --bdg -q 0.05 -f BAM
macs3 callpeak -t STAT1_30m_IFNa.bam -c INP_30m_IFNa_11E6.bam -n STAT1_30m_IFNa -g hs --bdg -q 0.05 -f BAM

# Obtain script to convert bedgraph to wiggle:
wget https://gist.githubusercontent.com/svigneau/8846527/raw/7fdf60b379245904f1aeb02427e098c06aeb670e/%2520bedgraph_to_wig.pl

# Change permissions to run scripts: 
chmod 755 bedgraph_to_wig.pl

# Convert bedgraph to wiggle:
perl bedgraph_to_wig.pl --bedgraph STAT1_6h_IFNa_treat_pileup.bdg --wig STAT1_6h_IFNa.wig --step 50
perl bedgraph_to_wig.pl --bedgraph STAT1_6h_IFNa_control_lambda.bdg --wig INP_6h_IFNa.wig --step 50

# Change permissions so scripts are executable:
chmod 755 bedGraphToBigWig
chmod 755 fetchChromSizes

# Fetch hg19 chrom sizes reference file:
./fetchChromSizes hg19 > hg19.chrom.sizes

# Convert bedgraphs to bigwigs:
./bedGraphToBigWig STAT1_6h_IFNa_treat_pileup.bdg hg19.chrom.sizes STAT1_6h_IFNa.bigwig
./bedGraphToBigWig STAT1_6h_IFNa_control_lambda.bdg hg19.chrom.sizes INP_6h_IFNa.bigwig

# Create alternate peaks.bed file without coumn 5:
cut -f1,2,3,4 STAT1_6h_IFNa_peaks.narrowPeak > STAT1_6h_IFNa_peaks.bed

less STAT1_6h_IFNa_peaks.bed
less STAT1_6h_IFNa_peaks.narrowPeak 
less STAT1_6h_IFNa_peaks.xls

# Load homer and make tag directories:
# makeTagDirectory STAT1_30m_IFNa_tagdir norm_STAT1_30m_IFNa.sam
# makeTagDirectory STAT1_6h_IFNa_tagdir norm_STAT1_6h_IFNa.sam
makeTagDirectory STAT1_30m_IFNa_tagdir STAT1_30m_IFNa.bam
makeTagDirectory STAT1_6h_IFNa_tagdir STAT1_6h_IFNa.bam

# Find differential peaks between 30min and 6h IFNα STAT1 ChIPs:
getDifferentialPeaks STAT1_30m_IFNa_peaks.narrowPeak STAT1_30m_IFNa_tagdir/ STAT1_6h_IFNa_tagdir/ > IFNa30m_diff_peaks.txt
getDifferentialPeaks STAT1_6h_IFNa_peaks.narrowPeak STAT1_6h_IFNa_tagdir/ STAT1_30m_IFNa_tagdir/ > IFNa6h_diff_peaks.txt

less IFNa6h_diff_peaks.txt

# Run CEAS gene-centered annotation of peaks with STAT1 6h IFNα peaks bed file:
ceas --name=STAT1_6h_IFNa_bed_CEAS -g hg19.refGene -b STAT1_6h_IFNa_peaks.bed

# Run CEAS gene-centered annotation of peaks with STAT1 6h IFNα peaks bed file:
ceas --name=STAT1_6h_IFNa_ISGs_CEAS --pf-res=50 -g hg19.refGene --dump --gn-groups=ISGs.txt --gn-group-names='ISGs' -w STAT1_6h_IFNa.wig

# Run sitepro analysis at STAT2 peaks for 6h IFNa STAT1 ChIP and Input:
sitepro --name=stat2peaks --pf-res=50 --span=3000 -b stat2_peaks.bed -w STAT1_6h_IFNa.wig -w INP_6h_IFNa.wig -l stat1_6h -l input_6h

# Run bedtools intersect on STAT1 and STAT2 peak locations:
bedtools intersect -wa -u -a STAT1_6h_IFNa_peaks.bed -b stat2_peaks.bed > STAT1_peaks_overlapping_STAT2.bed

# Count the lines in the output:
wc -l STAT1_peaks_overlapping_STAT2.bed

# Run findMotifsGenome.pl to search for de novo and known transcription factor motifs:
findMotifsGenome.pl STAT1_6h_IFNa_summits.bed hg19 STAT1_6h_IFNa_homer -size 200 -preparsedDir .

# http://great.stanford.edu/public/html/

# Compress Homer output folder:
tar -zcvf STAT1_6h_IFNa_homer.tar.gz STAT1_6h_IFNa_homer/

# Extract Homer folder:
tar -zxvf STAT1_6h_IFNa_homer.tar.gz

annotatePeaks.pl STAT1_6h_IFNa_peaks.bed hg19 > STAT1_6h_IFNa_peaks_annotatepeaks.txt -go STAT1_6h_IFNa_peaks_GO

less -S STAT1_6h_IFNa_peaks_annotatepeaks.txt 
