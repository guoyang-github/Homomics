# coding=utf-8
# guoyang
import subprocess
import argparse
import os.path
import sys
import os

usage='''
PyClone is a Bayesian clustering method for grouping sets of deeply sequenced somatic mutations 
into putative clonal clusters while estimating their cellular prevalences and accounting 
for allelic imbalances introduced by segmental copy-number changes and normal-cell contamination.

analysis of deeply sequenced (coverage > 100×) mutations to identify and quantify clonal populations 
in tumors, which extends to modeling mutations measured in multiple samples from the same patient.

The allelic prevalence of a mutation is a compound measure of several factors: 
the proportion of ‘contaminating’ normal cells, 
the proportion of tumor cells harboring the mutation and the number of allelic copies of the mutation 
in each cell, plus uncharacterized sources of technical noise.

four novel modeling advances:
	First, it uses beta-binomial emission densities, which models data sets with more variance 
		in allelic prevalence measurements more effectively than a binomial model.
	Second, flexible prior probability estimates (‘priors’) of possible mutational genotypes are used, 
		reflecting how allelic prevalence measurements are deterministically linked to zygosity 
		and coincident copy-number variation events.
	Third, Bayesian nonparametric clustering is used to discover groupings of mutations 
		and the number of groups simultaneously. This obviates fixing the number of groups a priori 
		and allows for cellular prevalence estimates to reflect uncertainty in this parameter.
	Fourth, multiple samples from the same cancer may be analyzed jointly to leverage the scenario 
		in which clonal populations are shared across samples.

http://pyyaml.org/wiki/PyYAMLDocumentation
###############################################
Usage
https://bitbucket.org/aroth85/pyclone/wiki/Usage


Overview

To run a PyClone analysis you need to perform several steps.
    1. Prepare mutations input file(s).
        1. Prepare .tsv input file.
        2. Run PyClone build_mutations_file --in_files TSV_FILE where TSV_FILE is the input file you have created.
    2. Prepare a configuration file for the analysis.
    3. Run the PyClone analysis using the PyClone run_analysis --config_file CONFIG_FILE command. Where CONFIG_FILE is the file you created in step 2.
    4. (Optional) Plot results using the plot_clusters and plot_loci commands.
    5. (Optional) Build summary tables using the build_table command.

###############################################
Tutorial
https://bitbucket.org/aroth85/pyclone/wiki/Tutorial

输入：
Before you can use PyClone you will need the following information for each mutation in the sample.
    1. Allelic count data from a sequencing experiment. PyClone requires that you specify the number of reads overlapping the mutation which contain the reference allele and the number of reads which contain the variant allele.
    2. The copy number of the genomic region containing the mutation. PyClone can work with either predictions of total copy number or parental (allele specific) copy number. In general performance will be better if you can specify parental copy number.
	3. tumour content

ps:We cannot usually tell what allele is maternal or paternal, so we use minor and major copy number instead. The convention is that the major copy number is the larger of the two values.
ps:If you only have total copy number for the tumour, not parental copy number, you can set the minor_cn to 0 and the major_cn to the predicted total copy number. When we use the PyClone build_mutations_file file command you will need to set the flag --var_prior total_copy_number in this case. By default the command assumes parental copy number information is being passed.
ps:Tools which can predict parental copy number are ideal, and even better are tools which also provide an estimate of tumour content.

应用范围：
In principle whole genomes shotgun sequencing (WGSS) or exome capture sequencing data could be used. 
In practice the depth of these approaches will be to low for an accurate PyClone analysis. 
The preferred approach is to use deep sequence data acquired by targeted amplicon sequencing 
or custom capture arrays.  建议大于100X



问题：
这是一个例子：    https://bitbucket.org/aroth85/pyclone/wiki/Tutorial
Positions in which only one of the cases has a variant genotype (AB or BB) are included in this dataset. 
Conceptually this is equivalent to a sample with four populations of diploid clones, 
which share no mutations. Because we excluded sex chromosomes, 
the total copy number of all positions is 2. 
For this dataset we still need to get the parental copy number for mutation. 
This can be done since we have the predicted genotype of the variants (AB or BB). 
In the case the mutation is AB the major and minor copy numbers would both be 1. 
In the case the genotype was BB the major copy number would be 2 and the minor copy number would be one. （应该是0）
Below is an example of the first two rows of one of the input files, SRR385938.tsv. 
The input files are located under the tsv/ folder.

mutation_id 	ref_counts 	var_counts 	normal_cn 	minor_cn 	major_cn 	variant_case 	variant_freq 	genotype
NA12156:BB:chr2:175263063 	3812 	14 	2 	0 	2 	NA12156 	0.0036591741 	BB
例子里面假设几个亚克隆结构，通过单个亚克隆得到了基因分型，判断了allelic的拷贝数，比如BB -》0 2
但如果从整体上得到AB的基因型，某个提供该变异的亚克隆同样可以是BB型，这样得到的allelicCN是不可确定的
还不如干脆就用total的
猜想可以通过抽样／逼近被估计的几个变量： 1. 浓度  2. 亚克隆比例  3. allelic情况
'''

class innerSoftware:
	def __init__(self):
		self.python=''
		self.pyclone=''
		

class PYCLONE:
	def __init__(self, tumorIdList, cnvList, vcfList, outDir, software=innerSoftware()):
		self.tumorIdList = tumorIdList.split(',')
		self.num = len(self.tumorIdList)
		self.cnvList = cnvList.split(',')
		self.vcfList = vcfList.split(',')
		self.outDir = os.path.abspath(outDir)
		self.software=software
		sys.exit('Unmatched sample number!') if len(self.cnvList) != self.num or len(self.vcfList) != self.num else None

		
	def pyClone(self):
		code=[self.software.java1_8 + ' -jar ' + self.software.picard + ' BuildBamIndex',
			'INPUT=' + self.infile
			#'OUTPUT=File'  #The BAM index file. Defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai.
		]		
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		subprocess.check_call(self.pyClone(), shell=True)
		
	def printJob(self):
		return self.pyClone() + '\n\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument('--tumorIdList', help = "tumor id list, seperated by ','", required = True)
	parser.add_argument('--cnvList', help = "controlfreeC _CNV list, seperated by ','", required = True)
	parser.add_argument('--vcfList', help = "mutect vcf list, seperated by ','", required = True)
	parser.add_argument('--outDir', help = "work directory", required = True)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	pyclone=PYCLONE(args.infile)
	pyclone()