# ddRAD_NWM

# Hi Amely,
# I have created this Github repository so that we are able to keep track of everything we will do for our IPS project.

# This is the one I wil modify (temporary version). Unless I add it to the Master (Which is the final version) this will remain as temp.

# THIS IS THE MOST UP TO DATE PIPELINE THAT I AM WORKING ON. 
# I AM CURRENTLY TRYING BBDuk.
# IF YOU WANT WE CAN CREATE ANOTHER FILE WITH RESULTS


#ddRAD PIPELINE - Data management and analysis

###################################################
############## Logging to the TACC ################
###################################################
ssh ssh lmv498@ls5.tacc.utexas.edu

###################################################
###### Downloading the data from GSAF #############
###################################################

# To download data from GSAF through Lonestar you must incorporate the  BioITeam start-up script
# into your profile and then run the command "gsaf_download.sh" to download your data. To

# To set up the BioITeam environment on LS5, which gives you access to the $BI directory, you need to
# add a line to your .bashrc file. Run the following line on Lonestar 5:

/work/01863/benni/bioiteam_ls5_env.sh

# Then logout and log back in.

# Create the folder where you will download the data. Your data is downloaded in $SCRATCH, however
# you only keep files here if you are using them. Otherwise keep your data in the $WORK directory

mkdir folder_name gsaf_download.sh
"http://gsaf.s3.amazonaws.com/JA15191.SA15082.html?AWSAccessKeyId=AKIAJPRTJ4ZGKJIUA6XA&Expires=
1440955025&Signature=k23lCfRK8q351MeIuLsH9zw007c%3D"

# Data is download in FASTQ format (typical NGS format).  FASTQ are like regular Fasta file but it is
# a format for storing sequence as well as quality. FastQC files  come compressed so you should leave
# them like that as they save space and most tools handle .gzip. However If you want to unzip the
# data to look at it you will need to unzip them.

gunzip JA15191Pool_CCGCGT_L006_R1_001.fastq.gz gunzip JA15191Pool_CCGCGT_L006_R2_001.fastq.gz

# GSAF gives you paired end sequencing data in two matching fastq format files, contining reads for
# each end sequenced -- for example:
	# Sample_ABC_L005_R1.fastq and Sample_ABC_L005_R2.fastq. 
# To look at them:

head -4 JA15191Pool_CCGCGT_L006_R1_001.fastq

# The data looks like this:

@HWI-D00289:185:C6FE5ANXX:6:1101:1429:2158 1:N:0:CCGCGT
GCATGCATGCCACTGTAAGCACAGTGTTTTTTATAGTCACAGTGGAGGACTGATGTTGCATGTCCATTGAGTTTACTTTCAGGAAACTTTGTCCCTTAAG
AACATATCCCAGGAAACCTGTTTTTC +
<BB<BFFFFFBF/FFFFB/FFFBFFFFBBFB<FFFFBFFB/FFFFFB<FBFFFFFFFF//<FFFFFF/FFF//B///FB<F//</<<FF/FFB//<<//
BFFFBFBFF/<</7/FB<7<B/B/<//

# Line 1: read identifier, which describes the machine, flowcell, cluster and grid coordinate.
# Except for the barcode information, read identifiers will be identical for corresponding entries
# in the R1 and R2 fastq files. A space separates the name from extra read information: 
#	- End number (1 for R1, 2 for R2) - Two qualtiy fields (N = not QC failed)
# 	- Barcode
# Line 2:   sequence reported by the machine. Line 3: always '+' from GSAF (it can optionally include
# a sequence description) Line 4: string of Ascii-encoded base quality scores, one character per base
# in the sequence. For each base, an integer quality score = -10 log(probabilty base is wrong) is
# calculated, then added to 33 to make a number in the Ascii #printable character range.

# R1 starts with the barcode (5bp) and then is followed by the RE1 recongnition site (4bp). 
# R2 Starts with the RE2 recongnition site.

# Counting the number of reads

grep -c '^@' JA15191Pool_CCGCGT_L006_R1_001.fastq

# If you want to look the files, don't use the 'cat' command alone, because the files are giant.
# Instead, we can use the 'less' which shows only what can be seen int the window.

less JA15192Pool_TGCAAA_L006_R2_001.fastq


###################################################
########@######## Running FASTQC  #################
###################################################

module load fastqc
fastqc filename.fastqc.zip


### Transfer data from TACC to your computer
scp -r
lmv498@ls5.tacc.utexas.edu:/scratch/02202/lmv498/ddRAD/JA16240/SA16089/Pool1_S5_L001_R1_001_fastqc.
zip
/Users/linavalencia/Documents/PhD_Research/Methods/NGS/ddRAD_Valencia&Martins/LibraryPrep/SA16089

###################################################
################# Quality Filter###################
###################################################

# Once you check your FastQC results, specially per base sequence quality report, you decide whether
# you trimm or not. There are many pipelines out there for trimming either by base or for for
# overall quality. The important thing is that you pick a pipeline that deals with pair end data.
# After filtering the R1 and R2 based on quality, there are potentially less reads which suggests
# that some of the resulting reads will not have a matching R1 or R2 read. Because most of the 
# pipelines used for the downstream analysis required paired data, Fastx TOOLKIT doesnt work.

#### Nonetheless, using FastX_Toolkit:

module load fastx_toolkit 

# Only module it once. FASTX is a tool for processing FASTA/FASTAQ data files.
# Paste the following and save the commands file:

cat JA15191Pool_CCGCGT_L006_R1_001_adaptrimm.fastq | fastq_quality_filter -Q33 -q 20 -p 90 >
JA15191Pool_CCGCGT_L006_R1_001.trim cat JA15191Pool_CCGCGT_L006_R2_001_adaptrimm.fastq |
fastq_quality_filter -Q33 -q 20 -p 90 > JA15191Pool_CCGCGT_L006_R2_001.trim

# As GSAF results are 33-based we have to add the -Q33 The options -q 20 -p 90 mean that 90% or more
# of all bases within the read should have PHRED quality of at least 20 (i.e., probability of error
# 1% or less, 1/100 errors) PHRED quality=10*(-lg(P.error))

# Compare how many reads remained after trimming AND filtering?
grep @HWI JA15191Pool_CCGCGT_L006_R1_001_adaptrimm.fastq | wc -l grep @HWI
JA15191Pool_CCGCGT_L006_R1_001.trim | wc -l

grep @HWI JA15191Pool_CCGCGT_L006_R2_001.adaptrimm | wc -l grep @HWI
JA15191Pool_CCGCGT_L006_R2_001.trim | wc -l


# WHen you assess quality you have to think of:
#### 1. Overall sequence quality
#### 2. Per base sequence quality
#### 3. Adapter contamination



##########################################################
#####DO WE HAVE ADAPTOR CONTAMINATION IN OUR SAMPLES? ####
###########################################################

# You can check this either by looking at your FastQC by looking at your overrepresented sequences
# or you can simply do a query.
#
### P5 PCR primer/flowcell capture site: 5’ AATGATACGGCGACCACCGAGA 3’
### Read1 primer site (also sequence for sequencing primer): 5’TCTACACTCTTTCCCTACACGACGCTCTTCCGATCT 3’
# (if only using one index)
### Read2 primer site (RC is  sequencing primer): 5’ AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC 3’
### P7 PCR primer/flowcell capture site: 5 ‘ ATCTCGTATGCCGTCTTCTGCTTG  3’

##### Looking for contamination in R2 (Reverse contaminant), you look at R1 primer site and P5 adapter
## contamination. Because R2 is the reverse complement of R1, you need to look for the reverse
## complement of R1 primer and P5 adapter
## R1 contamination (use reverse complement):
grep -c AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGA JA15191Pool_CCGCGT_L006_R2_001.fastq 

## P5 contamination (use reverse complement)
grep -c TCTCGGTGGTCGCCGTATCATT JA15191Pool_CCGCGT_L006_R2_001.fastq 

##### Looking for R1 contamination (Forward contaminant)you look for R2 primer and P7 adapter.
## It is expected to have sometimes more of P7 contaminants than of R2,
## because it is possible that simply the sequences was so small that not even the whole R2 was
## sequenced
#Looking for R2 primer contaminant. 
grep -c AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC JA15191Pool_CCGCGT_L006_R1_001.fastq #337

## Looking for P7 contaminations
grep -c ATCTCGTATGCCGTCTTCTGCTTG JA15191Pool_CCGCGT_L006_R1_001.fastq #117

####Checking how does sequences look:
cat JA15191Pool_CCGCGT_L006_R2_001.fastq  | grep AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCT  > adapter


##### The following pipeline is based on the suggestions made by the developer of bbmap - Brian Bushnell

# 1. Adapter trimming (using BBDuk) 
# 1.b If reads have an extra base at the end (like 2x151bp reads versus 2x150bp),it should be trimmed 
# here with the "ftm=5" flag.  That will occur before adapter-trimming.

# 2. Contaminant filtering for synthetic molecules and spike-ins such as Phix (using BBDuk)
# 2b. Quality-trimming and/or quality-filtering (Only recommended if you have very low-quality data or are doing something
# very sensitive to low-quality data, like calling very rare variants).

# 3. Paired-read merging.  Optional; mainly for assembly, clustering, or insert-size calculation (using BBMerge).
# 3b. RQCFilter runs BBMerge on all paired libraries for insert-size calculation, and uses the "cardinality" flag to simultaneously
# calculate the approximate number of unique kmers in the dataset, which can help estimate memory needs for assembly.


### If your data comes in different libraries/files you can concatenate the data after checking for quality. BE sure you always 
# concatenate in the same order and always keep R1 in one file and R2 in the other file.


cat Pool1_S5_L001_R1_001.fastq.gz Pool1_S5_L002_R1_001.fastq.gz > Pool1_S5_R1_001.fastq.gz

##########################################################
##################### ADAPTER REMOVAL ####################
#################### QUALITY TRIMMING ####################
##########################################################

###### 1. USING BBDuk from BBTOOLS
## ADAPTER TRIMMING:
# It is recommended to do first adapter removal and then quality filter to be able to maximize the matches between contaminants and 
# adapter sequence. If you do quality first you will loose data and thus make more difficult the matching.
# HOwever, it is recommended to do quality filter after merging reads. Once reads are merged, the quality assigned to the new sequence
# is the mean of R1 and R2, thus it is expected that the overall sequence quality increases and thus you loose less data.

# Be sure the whole command is in a line. If your data is in different libraries/files in your params file you can simply add as 
# many commands as libraries/files and in your launcher file specify a core per job.
# Create within your folder and adapters file

# BBDuk operates on kmers in one of 4 modes: right-trimming, left-trimming, masking, or filtering.  
# The default is filtering - any read matching a reference kmer will be discarded. In ktrim=r mode, once a reference kmer is
# matched in a read, that kmer and all the bases to the right will be trimmed, leaving only the bases to the left; 
# this is the normal mode for adapter trimming. 
# k = kmer length (at least the length of the barcode).
# mink = minimum distance of kmer (in those cases where the contaminants is shorter than k)
# hdist= "hamming distance" - number of mismatches
# hdist2 = haming distance for short kmers. Short kmers have a higher probability of attaching to random sequences and yield high 
# false positives. Therefore for short kmers is better to have a low hdist to minimize false positives
# tbo = also trim adapters based on pair overlap detection using BBMerge (which does not require known adapter sequences)
# tpe = trim both reads to the same length (in the event that an adapter kmer was only detected in one of them).
# qtrim = Trim read ends to remove bases with quality below trimq. Performed AFTER looking for kmers.
# trimq = Regions with average quality BELOW this will be trimmed.
# minlength = Reads shorter than this after trimming will be discarded.  Pairs will be discarded if both are shorter.
# minavgquality = Reads with average quality (after trimming) below  this will be discarded.

 
# If you want to do both and quality trimming use this code:
bbduk.sh in1=./raw/Pool1_S5_L002_R1_001.fastq.gz in2=./raw/Pool1_S5_L002_R2_001.fastq.gz out1=Pool1_S5_L002_R1_001_cl.fq 
out2=Pool1_S5_L002_R2_001_cl.fq outm=Pool1_S5_L002_R1_001_uncl.fq outm2=Pool1_S5_L002_R2_001_uncl.fq 
outs=Pool1_S5_L002_sing.fq stats=Pool1_S5_L002_stats.txt ref=adapters_LMV.fa 
ktrim=r k=23 mink=11 hdist=3 hdist2=1 rcomp=t tbo tpe minlen=30 minavgquality=20 minlen=20 qtrim=rl trimq=20

# JUst adapter trimming
bbduk.sh in1=./raw/concatenated_raw/Pool1_S5_R1_001.fq.gz in2=./raw/concatenated_raw/Pool1_S5_R2_001.fq.gz out1=Pool1_S5_L002_R1_001_cl.fq 
out2=Pool1_S5_L002_R2_001_cl.fq outm=Pool1_S5_L002_R1_001_uncl.fq outm2=Pool1_S5_L002_R2_001_uncl.fq 
outs=Pool1_S5_L002_sing.fq stats=Pool1_S5_L002_stats.txt ref=adapters_LMV.fa 
ktrim=r k=23 mink=11 hdist=3 hdist2=1 rcomp=t tbo tpe minlen=30 

##############################
###LOOKING AT OTHER OPTIONS###
##############################

#### 1. CUTADAPT:
# -m 20: discard any sequence that is smaller than 20 bases after trimming. This avoids problems trying to map very short, highly ambiguous sequences.
#  -O 10: only trim 3' adapter sequences that at least the first 10 bases of the adapter are seen. Prevents trimming short 3' sequences that just happen by chance to match the first few adapter sequence bases.
# -q 20: quality-trim the first read in each pair with a threshold of 20
# Allows three mistmaches between adapter and read.

# First trim the forward read R1, writing output to temporary files.
$BI/bin/cutadapt-1.3/bin/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --minimum-length 20 -O 10 -o JA15191Pool_CCGCGT_L006_R1_001temp.fastq -p JA15191Pool_CCGCGT_L006_R2_001temp.fastq ./raw/JA15191Pool_CCGCGT_L006_R1_001.fastq ./rawJA15191Pool_CCGCGT_L006_R2_001.fastq

#Then trim the reverse read R2, using the temporary files as input:
$BI/bin/cutadapt-1.3/bin/cutadapt  -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCT --minimum-length 20 -O 10 -o $SCRATCH/ddRAD/JA15191/JA15191Pool_CCGCGT_L006_R2_001_cutadapt.fastq -p $SCRATCH/ddRAD/JA15191/JA15191Pool_CCGCGT_L006_R1_001_cutadapt.fastq $SCRATCH/ddRAD/JA15191/JA15191Pool_CCGCGT_L006_R2_001temp.fastq $SCRATCH/ddRAD/JA15191/JA15191Pool_CCGCGT_L006_R1_001temp.fastq

#Finally, remove the temporary files:
rm JA15191Pool_CCGCGT_L006_R2_001temp.fastq JA15191Pool_CCGCGT_L006_R2_001temp.fastq

### With the new version of cutadapt not available in the TACC you do both reads at the same time
#cutadapt - -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA - A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCT  -o $SCRATCH/ddRAD/JA15191/JA15191Pool_CCGCGT_L006_R1_001_adaptrimm.fastq 
#$SCRATCH/ddRAD/JA15191/JA15191Pool_CCGCGT_L006_R2_001_adaptrimm.fastq $SCRATCH/ddRAD/JA15191/JA15191Pool_CCGCGT_L006_R1_001.fastq $SCRATCH/ddRAD/JA15191/#JA15191Pool_CCGCGT_L006_R1_001.fastq

#### 2. TRIMMOMATIC:
# Verify that the adapter file has your desired adapters. Be sure that you have both prefix 1 and 2,
# read primers and flowcell primer.
### Illumina clip: Trim adapters. Uses the adapter file you specify.
# Trimmomatic will look for seed matches (16 bases) allowing maximally 3 mismatches. These seeds
# will be extended and clipped. 
### For palindrome clipping If a score of 30 (threshold) is reached, about 50 bases must matched between
# the two adapter ligated reads (prefix), or in the case of single ended reads (simple clip threshold) a score of 10,
# about 18 bases must match between read and adapter.
# Remove low quality or N bases (below quality 20) at the beggining (TRAILING) and end (HEADING) of
# the read.
# Scan the read with a 5-base wide sliding window, cutting when the average quality per base drops
# below 20
# Drop reads which are less than 30 bases long after these steps 
# 2 input files and 4 output files (2 for matched R1 and R2 and 2 for unmatched R1 and R2)
 
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip 

nano
trimming

#copy the info below and be sure is only on one line in TACC.
java -jar /home1/02202/lmv498/bin/Trimmomatic-0.33/trimmomatic-0.33.jar PE -phred33 -trimlog
trimlog_Pool1_S5_L001.txt ../raw/Pool1_S5_L001_R1_001.fastq ../raw/Pool1_S5_L001_R2_001.fastq
Pool1_S5_L001_R1_paired.fastq Pool1_S5_L001_R1_unpaired.fastq
Pool1_S5_L001_R2_paired.fastq Pool1_S5_L001_R2_unpaired_TEST.fastq
ILLUMINACLIP:/home1/02202/lmv498/bin/Trimmomatic-0.33/adapters/TruSeq3-PE-2_LMV.fa:3:30:10
LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:30





##########################################################
############### CONTAMINANT FILTERING ####################
#################### PHIX VIRUS ########3#################
##########################################################

# The PhiX genome is used in library prepartion and is used as a control during Illumina sequencing.
# It provides a quality control for cluster generation, sequencing, and alignment, and a calibration control for 
# cross-talk matrix generation, phasing, and prephasing.
# It is a balanced and diverse library that can help mitigate sequencing challenges in unbalanced and low diversity libraries.
# It serves as a control for unbalanced samples (AT or GC content less than 40%), for low diversity samples, and for troubleshooting
# cluster generation problems (allows to determine if an error is related to library preparation or not).
# Low diversity libraries are libraries where a significant number of the reads have the same sequence. The base composition 
# shifts because the reads are no longer random.
# The adapter is the genome of the virus

#To remove phix:
bbduk.sh in1=./trimmed_adapteronly/Pool1_S5_R1_001_cl.fq in2=./trimmed_adapteronly/Pool1_S5_R2_001_cl.fq out1=Pool1_S5_L002_R1_001_clphix.fq 
out2=Pool1_S5_L002_R2_001_clphix.fq outm=Pool1_S5_L002_R1_001_unclphix.fq outm2=Pool1_S5_L002_R2_001_unclphix.fq 
outs=Pool1_S5_L002_singphix.fq stats=Pool1_S5_L002_statsphix.txt ref=/home1/02202/lmv498/bin/bbmap/resources/phix174_ill.ref.fa 
k=31 hdist=1 stats=stats.txt



##########################################################
############### DEMULTIPLEXING ###########################
##########################################################
### USING deML
# Computes the likelihood of assignment of a read to all potential samples of origin, assign each read to the most likely sample
# and compute the uncertainty of the assignment. It identifies all possible index sequences (barcodes) from a user provided list (ID)
# within a given number of mismatches, and estimates the likelihood of having sequenced the index given that it comes from a single sample.

# To be able to run deML you need three files:
# 1. Index file: First 5bp (assumed barcode) of your R1
# 2. R1 without the first 5bp
# 3. R2 

# Use FASTX tools to trim the first 5bp (resultind in input R1) and the rest of the read (resulting in INDEX)

# -f 6 : first base to keep
fastx_trimmer -Q33 -f 6  -i Pool1_S5_R1_001_clphix.fq -o ../trimmed/Pool1_S5_R1_001.fq
# -l 5 : Last base to keep
fastx_trimmer -Q33 -l 5  -i Pool1_S5_R1_001_clphix.fq -o ./trimmed/Pool1_S5_R1_001_ID.i

# However deML does not recognize the read name of GSAF output files:
# @K00179:39:H7JCJBBXX:1:1101:6837:1877 1:N:0:CGTACG
# And thus it needs to be change to the following:
# @K00179:39:H7JCJBBXX:1:1101:6837:1877/1

# To change names:
# For R1: 
sed -i.bak 's/\s.*/\/1/' Pool1_S5_R1_001_fnl.fq
# For R2:
sed -i.bak 's/\s.*/\/2/' Pool1_S5_R2_001_clphix.fq
# For ID:
sed -i.bak 's/\s.*//' Pool1_S5_R1_001_ID.fq

# Now you have three file ready to be runned in deML.

# Create a barcode.txt file

#Index1	Name
ATCGT	sample

#Run deML
deML -i barcodes.txt -f Pool1_S5_R1_001.fq -r Pool1_S5_R2_001.fq -if1 Pool1_S5_R1_001_ID.fq -o Pool1_S5_R1_001_fracconf.fq -s stats_fracconf.txt
deML -i barcodes.txt -f Pool1_S5_R1_001.fq -r Pool1_S5_R2_001.fq -if1 Pool1_S5_R1_001_ID.fq --mm 2 -o Pool1_S5_R1_001_fracconf.fq -s stats_fracconf.txt
deML -i barcodes.txt -f Pool1_S5_R1_001.fq -r Pool1_S5_R2_001.fq -if1 Pool1_S5_R1_001_ID.fq --fracconf 5 -o Pool1_S5_R1_001_fracconf.fq -s stats_fracconf.txt
deML -i barcodes.txt -f Pool1_S5_R1_001.fq -r Pool1_S5_R2_001.fq -if1 Pool1_S5_R1_001_ID.fq --fracconf 5 --mm 1 -o Pool1_S5_R1_001_fracconf.fq -s stats_fracconf.txt

## I run all parameters and using mm=2 and franccof=5 gives the greatest number of reads. I think mm=2 is too much so I decided to use the data
# that has only one mismatch using default parameters (not using franccof=5).

#Options: Related to the probability of assigning a sequenced Index to a particular barcode. If failed -- UNKNOWN
#--rgqual (Z0): Likelihood of having sequenced the index given that it comes from a given sample. 
#--fracconf (Z1): Compares the Z0 score for the best read group to all the other Z0 scores. Probability of missasignment (if a read equaly
# matches to barcodes. Make sure that there are two read groups with comparable Z0 to avoid misassignments. If failed -- CONFLICT
#--wrongness (Z2): Probability of mispairing (for double indexing). If failed -- WRONG
#--mm  Maximum # of tolerated mismatches [2] 

###deML 
# "unknown": More likely that you belong to that sample than any other but that probability is low, so you fail
# "conflict": More likely that you belong to that sample than any other but the probability that you belong to another sample 
# is almost equally likely, so you fail
# "wrong": It's more likely that you belong to that sample than any other but it seems your indices (in the case of double 
# indexing) are mispaired so you fail. NOT OUR CASE.
# "assigned": Congrats! You have none of the problems listed above

##TO RENAME FILE from PoolD_WENTAMO.fq_WENTAMO_r1.fq.gz TO  WENTAMO_r1.fq.gz
for file in *.fq_*; do echo mv $file ${file//*.fq_/ }; done
for file in *gz; do mv ${file//_r/_R}; done 

# 3. BBMerge: (instead of PEAR): Merges overlapped pairs.  This is a result of having very small reads, and consequently having
# R1 and R2 overlapping. In our case however, it is expected to have a great proportion of overlapped reads as we have 300bp reads and
# did 2X150.bp.
# WHY? - Single reads are easier to assemble, and longer reads give a better assembly;
# 	   - Merging can detect or correct errors where the reads overlap, yielding a lower error rate than the raw reads. 
# However, if the reads are merged incorrectly – yielding a single read with a different length than the actual insert size of
# the original molecule – errors will be added.

                     
bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>

### TO DO BBMERGE IN LOOP
INPUT_DIR=/scratch/02202/lmv498/ddRAD/JA16240/SA16089/quality_filter/BBDUK/demultiplexed/deML_mm1
OUTPUT_DIR=/scratch/02202/lmv498/ddRAD/JA16240/SA16089/quality_filter/BBDUK/merged

cd /scratch/02202/lmv498/ddRAD/JA16240/SA16089/quality_filter/BBDUK/demultiplexed/deML_mm1/
for file in *_r1.fq.gz; do echo "bbmerge.sh in1=$INPUT_DIR/$file in2=$INPUT_DIR/${file//r1/r2}  out=$OUTPUT_DIR/${file//_r1/_merged} outu1=$OUTPUT/${file/_r1/_r1_unmerged} 
outu2=$OUTPUT_DIR/${file/_r1/_r2_unmerged}" >> $OUTPUT_DIR/bbmerge.params; done 

# ALTERNATIVELY create a direct link to input files:
ln -s PATH fastq_mm1  
mkdir paired
mkdir merged
# Within you input (deML files) folder

for file in *R1*gz; do echo bbmerge.sh in=./fastq_mm1/$file in2=./fastq_mm1/${file//R1/R2} out=./merged/${file//R1/} outu1=./paired/$file outu2=./paired/%${file//R1/_R2} minoverlap=10 "&>"${file//_R1.f1.gz/.bbmerge.log) >> bbmerge_mm1.params ;done

## TO RUN PYRAD
#. 1 Create a params.txt specifying . We will create two .txt depending on the numbers of processors (## 9.) that will be used.
# params.txt will be used to run s3 as it is very time consuming
# params_n48. txt  for all steps.

#Step 2: Run all individuals at the same time in one node and one core.
# params_n48.txt, as you are running it in 48 processors.
## 9. NQual: max # sites with qual < 20 (or see line 20)(s2)  --> 10
## 20.opt.: phred Qscore offset (def= 33)  --> 33
## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict  --> 1

#Step 3:
# Run in paralel, running each individual at a time.
# Do a loop to create a params.txt for each of your samples (MERGED AND UNMERGED?), specifying the prefix in #15.
# To do so:
#1. Create a inds.list (with the samples names).

#2. Create params files

mkdir params
export INDS=$(cat inds.list)
echo $INDS
for ind in $INDS; do cp params.txt params/"params_"$ind".txt"; done

#open Excel file and replace samples names with your own samples.
#Copy the commands info into a .sh and run it
# You will have two .sh files, one to change the prefix and one to change the output.
bash command.sh
# Within the params folder, now you have a pyrad params file for each one of your samples.

# To create your .params file to run do:

cd params 
for file in *.txt; do echo pyrad -p "params/"$file -s 3 "&&"  pyrad -p "params/"$file -s 2 "&&"  pyrad -p "params/"$file -s 3 "&&"  pyrad -p "params/"$file -s 2 "&&"  pyrad -p "params/"$file -s 3 "&&"  pyrad -p "params/"$file -s 3 "&&"  pyrad -p "params/"$file -s 3 >> pyrad_s3.params;done #si lo hace una vez pasa a la siguiente muestra
mv pyrad_s3.params
# we have specified the command 10 times, beacase is some cases it doesnt read the line and skip the sample. So we want to avoid that

 # I am currently running step 2 and will soon run step 3 for all samples
 # I am also running our old pipeline to compare the results. I will send you a summary of results soon.
 

