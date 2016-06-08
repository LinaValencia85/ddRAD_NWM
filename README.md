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

##########################################################
##################### ADAPTER REMOVAL ####################
#################### QUALITY TRIMMING ####################
##########################################################

###### 1. USING BBDuk from BBTOOLS
## ADAPTER TRIMMING:
# I had a lot of issues running it. Be sure the whole command is in a line.
# Create within your folder and adapters file

bbduk.sh -Xmx1g in1=./raw/Pool1_S5_L002_R1_001.fastq.gz in2=./raw/Pool1_S5_L002_R2_001.fastq.gz out1=trimmed_Pool1_S5_L002_R1_001.fq out2=trimmed_Pool1_S5_L002_R2_001.fq outm=unclean_Pool1_S5_L002_R1_001.fq outm2=unclean_Pool1_S5_L002_R2_001.fq outs=singletons_Pool1_S5_L002.fq stats=stats_Pool1_S5_L002.txt ref=adapters_LMV.fa ktrim=r k=23 mink=11 hdist=1 rcomp=t minkmerfraction=0.5 tpe tbo minavgquality=20 maxns=-1 minlen=20 qtrim=rl trimq=20


# ktrim=fl  		Trim reads to the left remove bases matching reference kmers.
# k = 				kmer length (at least the length of the barcode).
# hdist= 			means "hamming distance" and refers to number of mismatches
#tbo				trim adapters based on pair overlap detection (which does not require known adapter sequences)
#tpe				Trim both reads to the same length (in the event that an adapter kmer was only detected in one of them).
#qtrim=rl   		Trim read ends to remove bases with quality below trimq. Performed AFTER looking for kmers.
#trimq=6    		Regions with average quality BELOW this will be trimmed.
#minlength=10   	Reads shorter than this after trimming will be discarded.  Pairs will be discarded if both are shorter.
# minavgquality=0   Reads with average quality (after trimming) below  this will be discarded.

Subsampling?
reads=-1            If positive, quit after processing X reads or pairs.
samplerate=1        Set lower to only process a fraction of input reads.
 

