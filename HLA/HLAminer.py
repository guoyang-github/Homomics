# coding=utf-8
import subprocess
import argparse
import os.path
import sys
import os
import shutil
# 20170829
#注意数据库升级



usage='''
Manual Reference Pages – HLAminer - Derivation of HLA class I and class II predictions from shotgun sequence datasets
* This manual assumes that you have a working knowledge of unix, and some shell and perl scripting experience

NAME
  HLAminer - Derivation of HLA class I and class II predictions from shotgun sequence datasets


应用
RNA-Seq or WGS and exon capture

HLA-I类分子：内源性抗原的递呈分子
HLA-Ⅱ类分子：外源性抗原的递呈分子

--------
SYNOPSIS 示例
========

  For RNAseq:
  1. Copy ./test-demo/    eg. cp -rf test-demo foo
  2. In folder "foo", edit the patient.fof file to point to your NGS RNAseq data.  Ensure all paths are ok.
  3. For HLA Predictions by Targeted Assembly of Shotgun Reads: execute ./HLAminer/foo/HPTASRrnaseq.sh 
     For HLA Predictions by Read Alignment: execute ./HLAminer/foo/HPRArnaseq.sh



OVERVIEW
========

Derivation of HLA class I and class II predictions from shotgun sequence datasets (HLAminer) by:
1) Targeted Assembly of Shotgun Reads (HPTASR)
2) Read Alignment (HPRA)
两种方法 1.拼接的方法 2.直接比对的方法

DESCRIPTION
===========

The HLA prediction by targeted assembly of short sequence reads (HPTASR), 
performs targeted de novo assembly of HLA NGS reads and align them to reference HLA alleles 
from the IMGT/HLA sequence repository using commodity hardware with standard specifications (<2GB RAM, 2GHz).  
Putative HLA types are inferred by mining and scoring the contig alignments and an expect value is determined for each.  
The method is accurate, simple and fast to execute and, for transcriptome data, requires low depth of coverage. 
Known HLA class I/class II reference sequences available from the IMGT/HLA public repository are read by TASR using default options (Warren and Holt 2011) 
to create a hash table of all possible 15 nt words (k-mers) from these reference sequences. 
Note that this parameter is customizable and larger k values will yield predictions with increased specificity (at the possible expense of sensitivity). 
Subsequently, NGS data sets are interrogated for the presence of one of these kmers (on either strand) at the 5’ or 3’ start. 
Whenever an HLA word is identified, the read is recruited as a candidate for de novo assembly. 
Upon de novo assembly of all recruited reads, a set of contigs is generated.  
Only sequence contigs equal or larger than 200nt in length are considered for further analysis, as longer contigs better resolve HLA allelic variants.  
Reciprocal BLASTN alignments are performed between the contigs and all HLA allelic reference sequences. 
HPTASR mines the alignments, scoring each possible HLA allele identified, 
computing and reporting an expect value (E-value) based on the chance of contigs characterizing given HLA alleles and, 
reciprocally, the chance of reference HLA alleles aligning best to certain assembled contig sequences

The HLA prediction from direct read alignment (HPRA) method is conceptually simpler and faster to execute, 
since paired reads are aligned up-front to reference HLA alleles.  
Alignments from the HPTASR and HPRA methods are processed by the same software (HLAminer.pl) 
to derive HLA-I predictions by scoring and evaluating the probability of each candidate bearing alignments.




操作步骤总结：
=======
1. 解压
2. bin文件夹下所有perl shell等脚本及内部程序bwa等的软件路径修正
3. 配置ncbiBlastConfig，软链接到所用的blast版本
4. 配置patient.fof，写入fq data路径
5. 运行shell脚本开始执行（如HPTASRwgs_classI.sh）



INSTALL安装及配置
#####################
1. Download and decompress the tar ball
gunzip HLAminer_v1-3.tar.gz 
tar -xvf HLAminer_v1-3.tar
2. Make sure you see the following directories:
./bin
./databases
./docs
./test-demo

3. Read the docs in the ./docs/ folder
4. Change/Add/Adjust the perl shebang line of each .pl and .sh script in the ./bin/ folder as needed

From direct Read Alignment (HPRA, faster but less accurate):
HPRArnaseq_classI.sh
HPRArnaseq_classI-II.sh
HPRAwgs_classI.sh
HPRAwgs_classI-II.sh
-and for single end reads-
HPRArnaseq_classI_SE.sh
HPRArnaseq_classI-II_SE.sh
HPRAwgs_classI_SE.sh
HPRAwgs_classI-II_SE.sh


From Targeted Assembly (HPTASR, longer but more accurate):
HPTASRrnaseq_classI.sh
HPTASRrnaseq_classI-II.sh
HPTASRwgs_classI.sh
HPTASRwgs_classI-II.sh

*Running HPTASRwgs(rnaseq)_classI-II.sh will take longer than HPTASRwgs(rnaseq)_classI.sh, due to the reciprocal BLAST step.  
You may remove this step from the former (and HLAminer.pl command) to speed things up.  However, this step is helpful in weeding out spurious alignments to HLA references.  
That said, if you're solely interested in HLA-I, you have the option to run the latter set of scripts [HPTASRwgs(rnaseq)_classI.sh].

Also, in the ncbiBlastConfig2-2-XX.txt files (bin and test-demo directories), you may adjust the number of threads and number of reported alignments to speed things up. 
The options have different name depending on the blast version, refer to the blast manual
eg.
v2.2.22
option:description
-a:threads
-v:number of descriptions
-b:number of alignments

v2.2.28
-num_threads:threads
-max_target_seqs:number of hit sequences to report (when output is 5/xml)

In our hands, a few tests show that blast 2.2.22 may be faster than blast+ (2.2.28) while
producing accurate results - HLAminer (Warren et al. 2012) was thoroughly tested with 2.2.22

NCBI blast may be downloaded from:
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/
-or-
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/


HLAminer.pl
parseXMLblast.pl
TASR

5. You must install perl module Bio::SearchIO to use HPTASR
6. Edit the fullpath location of bwa and other software dependencies in the shell scripts in the ./bin/ folder, as needed
7. For your convenience, ncbi blastall and formatdb have been placed in the ./bin/ folder and executed from the following shell scripts:

NAME,PROCESS,NGS DATA TYPE,PREDICTIONS
HPRArnaseq_classI.sh,Paired read alignment,RNAseq (transcriptome),HLA-I A,B,C genes
HPRArnaseq_classI-II.sh,Paired read alignment,RNAseq (transcriptome),HLA-I A,B,C and HLA-II DP,DQ,DR genes

HPRAwgs_classI.sh,Paired read alignment,Exon capture (exome) and WGS (genome),HLA-I A,B,C genes
HPRAwgs_classI-II.sh,Paired read alignment,Exon capture (exome) and WGS (genome),HLA-I A,B,C and HLA-II DP,DQ,DR genes

HPTASRrnaseq_classI.sh,Targeted assembly of sequence reads,RNAseq (transcriptome),HLA-I A,B,C genes
HPTASRrnaseq_classI-II.sh,Targeted assembly of sequence reads,RNAseq (transcriptome),HLA-I A,B,C and HLA-II DP,DQ,DR genes

HPTASRwgs_classI.sh,Targeted assembly of sequence reads,Exon capture (exome) and WGS (genome),HLA-I A,B,C genes
HPTASRwgs_classI-II.sh,Targeted assembly of sequence reads,Exon capture (exome) and WGS (genome),HLA-I A,B,C and HLA-II DP,DQ,DR genes

Run those scripts by specifying the relative path ../bin/blastall or ../bin/formatdb in the shell scripts AND config file "ncbiBlastConfig.txt".  
Make sure HPRA and HPTASR are running in same-level directories to ../bin/ (eg. ../test-demo/)

8. Before running on your data, inspect the ./test-demo/ folder and familiarize yourself with the files and the execution.  
9. When you are ready and the demo works well, place the fullpath location of your short read fastq or fasta files in the "patient.fof" file.
10. Make sure that the following files are in your working directory:

patient.fof
ncbiBlastConfig.txt (specific to the version of blast you are using, see in
../bin and ../test-demo directories



ADDITIONAL INSTALL NOTES ON MAC OSX

install bioperl
---------------
sudo perl -MCPAN -e shell
install CJFIELDS/BioPerl-1.6.923.tar.gz
change shebang line in all PERL (.pl) scripts and location of bioperl on your system in HLAminer.pl

install ncbi blast
------------------
download and install blast-2.2.22-universal-macosx.tar.gz
change path to blast in ncbiconfig.txt

install homebrew
----------------
ruby -e "$(curl -fsSL
https://raw.githubusercontent.com/Homebrew/install/master/install)"

Since databases were indexed with older
version, had to re-index:
bwa index -a is HLA_ABC_CDS.fasta

change path to bwa in HPRA* shell scripts


COMMANDS AND OPTIONS
====================
可以修改HLAminer.pl参数
The shell scripts are set to filter out short (<200) contigs that would blur HPTAR predictions.  Feel free to adjust as you see fit.

Likewise, HLAminer.pl runs with the set defaults:
-z minimum contig size.......................<200> (HPTASR)
-i minimum % sequence identity...............<99>  (HPTASR / HPRA)
-q minimum log10 (phred-like) expect value...<30>  (HPTASR / HPRA)
-s minimum score.............................<1000> (HPTASR / HPRA)
-n consider null alleles (1=yes/0=no)........<0> (HPTASR / HPRA)

The minimum sequence identity applies to the short read paired alignment or blast alignment, depending on the choice made.  HLA predictions with a phred-like expect value lower than -q or a score lower than -s will not be diplayed.  Because IMGT/HLA reports numerous null alleles, an option exist to consider or not these unexpressed alleles. 

Likewise, TASR-based (Targeted Assembly of Short Reads) predictions could be
improved by using larger k values for assembly (-k). Experimentation for
choosing the ideal k to use depends on the input read length and is warranted. 


DATABASES数据库升级
=========

Follow these instructions to download updated HLA sequences from ebi/imgt 
(shell scripts to automatically download and format the databases exist in ./database/) 
and refer to README.txt in the ./database directory:


1) Coding HLA sequences

HLA CDS sequences from:

wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/fasta/A_nuc.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/fasta/B_nuc.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/fasta/C_nuc.fasta
cat A_nuc.fasta B_nuc.fasta C_nuc.fasta | perl -ne 'chomp;if(/\>\S+\s+(\S+)/){print ">$1\n";}else{print "$_\n";}' > HLA_ABC_CDS.fasta
../bin/formatdb -p F -i HLA_ABC_CDS.fasta
/home/pubseq/BioSw/bwa/bwa-0.5.9/bwa index -a is HLA_ABC_CDS.fasta

2) HLA genomic sequences 

To make the HLA genomic sequence database, execute these unix commands:

wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/fasta/A_gen.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/fasta/B_gen.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/fasta/C_gen.fasta
cat A_gen.fasta B_gen.fasta C_gen.fasta | perl -ne 'chomp;if(/\>\S+\s+(\S+)/){print ">$1\n";}else{print "$_\n";}' > HLA_ABC_GEN.fasta
../bin/formatdb -p F -i HLA_ABC_GEN.fasta
/home/pubseq/BioSw/bwa/bwa-0.5.9/bwa index -a is HLA_ABC_GEN.fasta

FOR YOUR CONVENIENCE, A SINGLE SHELL SCRIPT CAN BE RUN FROM ../database TO
UPDATE ALL IMGT-HLA SEQUENCE DATABASES AND CREATE BWA/BLAST INDEXES. JUST KEEP
IN MIND THAT THIS IS DONE WITH A SPECIFIC VERSION OF THESE TOOLS, IF YOU
UPGRADE, YOU WILL HAVE TO REGENERATE THE INDEXES
******************************
***../database/updateAll.sh***
******************************


3) P designation files

HLA alleles having nucleotide sequences that encode the same protein sequence 
for the peptide binding domains (exon 2 and 3 for HLA class I and exon 2 only for HLA class II alleles) 
will be designated by an upper case ‘P’ which follows the 2 field allele designation of the lowest numbered allele in the group.


Upgrade the P designation info from:

info:
http://hla.alleles.org/wmda/index.html

file:
http://hla.alleles.org/wmda/hla_nom_p.txt


OUTPUT FILES   结果展示
============

HLA predictions from read pair alignments:

HLAminer_HPRA.log
HLAminer_HPRA.csv

HLA predictions from targeted assemblies:

HLAminer_HPTASR.log
HLAminer_HPTASR.csv

The .log file tracks the process of HLA mining. It contains the following information:
-HLAminer command and parameters utilized
-Contig/read pair alignment output and best HLA hit for each
-Initial gene summary, score and expect value
-Final summary, listing all predictions by highest score (more likely).

The .csv file contains HLAminer predictions.  Predictions are listed by HLA
gene and ranked by highest score.  Predictions 1) and 2) are expected to
represent the two alleles for each.

eg.
----------------------------------------------------------------------
SUMMARY
MOST LIKELY HLA-I ALLELES (Confidence (-10 * log10(Eval)) >= 30, Score >= 500)
Allele,Score,Expect (Eval) value,Confidence (-10 * log10(Eval))
----------------------------------------------------------------------

HLA-A
Prediction #1 - A*26
        A*26:33,4179,5.22e-124,1232.8

Prediction #2 - A*33
        A*33:24,1791,2.41e-75,746.2

Prediction #3 - A*68
        A*68:05,597,1.85e-10,97.3


From these predictions, the individual is expected to be heterozygous with HLA-I A alleles A*26 and
A*33. The chance, expect value and confidence are different representations of the same
metric. The Confidence represents the Expect value - Eval - as a score in a
manner analoguous to the phred score employed in sequencing to quickly assess
the likelihood of a base being correct. 

Predictions/read pair are ambiguous when there are multiple predicted allele groups and/or protein coding alleles with the same score.



What's new in version 1.3?
==========================

A more concise HLA allele summary in HLAminer_HPTASR.csv and HLAminer_HPRA.csv (associated .log is unchanged and lists all predictions)
Keeps top two [highest-scoring by HLA group] predictions per gene and only the 'P' designated allele when the summary include HLA Sequences reported to have the same antigen binding domain.
For the original output, refer to the HLAminer_v1-2.pl included in the ./bin directory
A prediction example from MCF-7 PacBio RNA-seq reads is also provided


What's new in version 1.2?
==========================

Updated all HLA sequence databases
Corrected shell script that download HLA sequences to reflect change of location at EBI (ie. fasta sub folder) 
Added support for predictions from direct alignment of single-end reads

'''


class innerSoftware:
	def __init__(self):
		self.hlaminer='/HLAminer_v1.3'
		self.python='python'
		

class HLAMINER:
	def __init__(self, infile, type, projectdir, software=innerSoftware()):
		self.infile=os.path.abspath(infile)
		self.type=type
		self.projectdir=os.path.abspath(projectdir)
		self.software=software
		
	def HLAMiner(self):
		shutil.copy(self.software.hlaminer, self.projectdir)
		self.rundir=os.path.join(self.projectdir, 'HLAminer_v1.3', 'HLAminerNovo')
		shutil.copy(self.infile, self.rundir)
		if self.type == 'RNAseq':
			code = self.software.python + ' ' + os.path.join(self.rundir, HLAminerPwd.py) + ' ' + self.rundir + ' ' + 'HPTASRrnaseq_classI.sh'
		if self.type == 'WGS' or self.type == 'WES':
			code = self.software.python + ' ' + os.path.join(self.rundir, HLAminerPwd.py) + ' ' + self.rundir + ' ' + 'HPTASRwgs_classI.sh'
		#加这层python是为了改变chdir --> HLAminerPwd.py
		return code
		
	def __call__(self):
		subprocess.check_call(self.HLAMiner(), shell=True)
		
	def printJob(self):
		return self.HLAMiner() + '\n\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infiles',help="a file contains fastq files paths for read1 and read2 by line   patient.fof",required=True)
	parser.add_argument('--type', help="RNAseq/WGS/WES", choices=['RNAseq','WGS','WES'], required=True)
	parser.add_argument('--projectdir', help="project directory", default='.')
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	hlaMiner=HLAMINER(args.infile, args.type, args.projectdir)
	hlaMiner()