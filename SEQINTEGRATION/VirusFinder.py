# coding=utf-8
# guoyang
#20170627完成 未测试

#应用领域: WGS WES RNA-seq  可以用单端read数据，但本脚步强制双端
#功能: Report viruses   Report novel contigs   Report virus co-infections Output integration sites   Output consensus virus sequence   Rare variants in virus genome

# 4个步骤： 
# 1. Read subtraction    
# 2. Virus detection    
# 3. Virus integration detection    
# 4. Viral mutation detection

#准备工作 制造reference
#Host reference genome
#   bowtie2-build hg19.fa hg19
#   makeblastdb -in hg19.fa -dbtype nucl -out hg19
#Virus database (DB)
#       RINS package (Bhaduri et al., 2012, http://khavarilab.stanford.edu/resources.html)
#       Genome Information Broker for Viruses (GIB-V, http://gib-v.genes.nig.ac.jp/)
#       NCBI viral gene annotation (ftp://ftp.ncbi.nih.gov/refseq/release/viral/)
#   makeblastdb -in virus.fa -dbtype nucl -out virus

#流程：第一步抽出unmapped reads，第二步病毒检测可以有几种方法，1.给个病毒数据库里面有的序列号，直接用病毒库里的序列  2.给个已知的病毒序列，更直接
#3. 什么也不给，直接用unmapped的reads拼出几段contigs，再和病毒库去比较，看看最配哪条，就用病毒库里面那条完整序列当病毒reference（results-virus-top1.fa） 
# 4. 什么也不给，用unmapped的reads拼出几段contigs，也比不上病毒库，就成了novel contigs  
# 整个流程一次只能处理1个病毒，如果是co-infection的话，可以从第二步开始把病毒序列给换了，再做第三步第四步

#最终需要的输出文件： ‘virus.txt’, ‘virus-list.txt’, ‘contig.txt’, ‘integration-sites.txt’, and ‘viral-mutation.vcf’


#重要结果文件示例
#Chrommosome 1	Position 1	Strand 1	Chrommosome 2	Position 2	Strand 2	#Support reads (pair+softclip)	Confidence
#chr1	24,020,709	+	chrVirus	1	+	13+15	high

#The final column reports the confidence of VirusFinder in the actual position of virus-host integration. 
# High confidence means there are sufficient soft-clipped reads to support the virus integration locus. 
# Low confidence, however, indicates lack of soft-clipped reads for the accurate characterization of the locus.

import subprocess
import argparse
import os.path
import sys
import os

usage='''
Usage: VirusFinder.pl -c <configuration file> [options]
'''

configure='''
###########################################################################
#
# Configuration file for VirusFinder
#
###########################################################################


################################
## Input data can be: (a) an alignment file (in BAM format); or (b) FASTQ file(s) (preferred) – for 
## single end data, “fastq1” is required; for paired-end, both “fastq1” and “fastq2” are needed.
################################

# alignment_file = 
fastq1        = {fq1}
fastq2        = {fq2}

#mailto       = guoyang@novogene.com
thread_no     = {t}   #the number of threads for parallel computing

detect_integration = yes   # if no is provided, VirusFinder will not detect virus integrations
detect_mutation    = yes   # if no is provided, VirusFinder will not detect viral mutations


################################
## The full paths to the following third-party tools are required by VirusFinder:
################################

blastn_bin      = {software.blastn}
bowtie_bin      = {software.bowtie2}
bwa_bin         = {software.bwa}
trinity_script  = {software.trinity}
SVDetect_dir    = {software.SVDetect}


################################
## Reference files indexed for Bowtie2 and BLAST
################################

virus_database     = {database.virus}
bowtie_index_human = {database.bowtie_index_human}
blastn_index_human = {database.blastn_index_human}
blastn_index_virus = {database.blastn_index_virus}


##########################################
## Parameters of virus integration detection(VERSE algorithm). They are ignored for single-end data
##########################################

detection_mode     = sensitive   #Possible values: normal, sensitive; default value: normal. If not specified, VirusFinder runs in normal detection mode. 
flank_region_size  = 4000   ##Suggested values: >2000; default: 4000; if detection_mode = normal, it (and ‘sensitivity_level’ below) will be ignored.
sensitivity_level  = 1   ##Suggested values: 1~6; default value: 1; greater value means higher sensitivity, and accordingly more computation time.

##########################################
## Parameters of virus detection. Smaller "min_contig_length", higher sensitivity
##########################################

min_contig_length  = 300
blastn_evalue_thrd = 0.05
similarity_thrd    = 0.8 
chop_read_length   = 25
minIdentity        = 80
'''


class innerSoftware:
    def __init__(self):
        self.virusfinder = ''
        self.blastn = ''
        self.bowtie2 = ''
        self.bwa = ''
        self.trinity = ''
        self.SVDetect = ''
		
class innerDatabase:
    def __init__(self):
        self.virus = ''
        self.bowtie_index_human = ''
        self.blastn_index_human = ''
        self.blastn_index_virus = ''


class VIRUSFINDER:
    def __init__(self, fastq1, fastq2, virus, output, markdup, threads, software = innerSoftware(), database = innerDatabase()):
        self.fastq1 = os.path.abspath(fastq1)
        self.fastq2 = os.path.abspath(fastq2)
        self.virus = virus
        self.output = os.path.abspath(output)
        self.markdup = markdup
        self.threads = threads
        self.software = software
        self.database = database
		
    def VirusFinder(self):
        self.configureFile = os.path.join(self.output, 'virusfinder.config')
        with open(self.configureFile, 'w') as conf:
            conf.write(configure.format(database = self.database, software = self.software, fq1 = self.fastq1, fq2 = self.fastq2, t = self.threads))
        code=[self.software.virusfinder,
            '--config ' + self.configureFile,
            '--output ' + self.output,
            '--markdup ' + self.markdup
        ]
        if self.virus:
            code += ['--virus ' + self.virus]
        return '   \\\n\t'.join(code)
		
		
    def printJob(self):
        return self.VirusFinder() + '\n\n'

    def __call__(self):
        subprocess.check_call(self.VirusFinder(), shell=True)
		


def paramsParse():
    parser=argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('--fastq1', help = "read 1 fastq", required = True)
    parser.add_argument('--fastq2', help = "read 2 fastq", required = True)
    parser.add_argument('--virus', help = "The sequence file (in fasta format) of the virus. Not required", default = None)
    parser.add_argument('--output', help = "The directory to store software output, default is current working directory", default = '.')
    parser.add_argument('--markdup', help = "Mark duplicate reads. Accepted inputs: y[es], n[o]. Default value is y[es]. Duplicate reads will not be used for variant calling. For ultra-deep amplicon sequencing data, 'no' should be used for this argument.", default = 'y')
    parser.add_argument('--threads', help = "the number of threads for parallel computing", default = 8)
    return parser.parse_args()

if __name__ == '__main__':
    args=paramsParse()
    virusFinder=VIRUSFINDER(args.fastq1, args.fastq2, args.virus, args.output, args.markdup, args.threads)
    virusFinder()