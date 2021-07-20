# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar picard.jar MergeVcfs \
      I=input_1.vcf \
      I=input_2.vcf \
      O=merged_files.vcf
'''

class innerSoftware:
	def __init__(self):
		self.java1_8='java'
		self.picard='picard.jar'
		self.gatk='GenomeAnalysisTK.jar'
		self.bcftools=''
		
		
class innerDatabase:
	def __init__(self):
		self.reference='human_g1k_v37_decoy.fasta'

class MERGEVCF:
	def __init__(self, infile, outfile, type, software=innerSoftware(), database=innerDatabase()):
		self.infile=infile.strip().split(',')
		self.outfile=os.path.abspath(outfile)
		self.type=type
		self.software=software
		self.database=database
		
	def verticalMergeVCFpicard(self):
		'''
		Merges multiple VCF or BCF files into one VCF file. Input files must be sorted by their contigs and, within contigs, by start position. The input files must have the same sample and contig lists. 
		'''
		code=[self.software.java1_8 + ' -jar ' + self.software.picard + ' MergeVcfs',
			'OUTPUT=' + self.outfile
			#'SEQUENCE_DICTIONARY'  #The index sequence dictionary to use instead of the sequence dictionary in the input file Default value: null. 
		]
		code+=map(lambda x: 'INPUT=' + x, self.infile)
		return '   \\\n\t'.join(code)
	
	def verticalMergeVCFgatk(self):
		'''
		All the input VCFs (or BCFs) contain the same samples in the same order.
		The variants in each input file are from non-overlapping (scattered) intervals.
		'''
		code=[self.software.java1_8 + ' -cp ' + self.software.gatk + ' org.broadinstitute.gatk.tools.CatVariants',
			'-R ' + self.database.reference,
			'-out ' + self.outfile,
			'-assumeSorted'  #assumeSorted should be true if the input files are already sorted (based on the position of the variants)
		]
		code+=map(lambda x: '-V ' + x, self.infile)
		return '   \\\n\t'.join(code)	
		
	def horizontalMergeVCFgatk(self):
		'''
		CombineVariants reads in variants records from separate ROD (Reference-Ordered Data) sources and combines them into a single VCF. Any number of sources can be input. This tool aims to fulfill two main possible use cases, reflected by the two combination options (MERGE and UNION), for merging records at the variant level (the first 8 fields of the VCF) or at the genotype level (the rest).

		MERGE: combines multiple variant records present at the same site in the different input sources into a single variant record in the output. If sample names overlap, then they are "uniquified" by default, which means a suffix is appended to make them unique. Note that in version 3.3, the automatic uniquifying was disabled (unintentionally), and required setting `-genotypeMergeOptions UNIQUIFY` manually.
		UNION: assumes that each ROD source represents the same set of samples (although this is not enforced). It uses the priority list (if provided) to emit a single record instance at every position represented in the input RODs.
		'''
		code=[self.software.java1_8 + ' -jar ' + self.software.gatk,
			'-T CombineVariants',
			'-R ' + self.database.reference,
			'-o ' + self.outfile,
			'-genotypeMergeOptions UNIQUIFY'
			#UNIQUIFY
			#Make all sample genotypes unique by file. Each sample shared across RODs gets named sample.ROD.
			#PRIORITIZE
			#Take genotypes in priority order (see the priority argument).  must with the par '-priority',Ordered list specifying priority for merging
			#UNSORTED
			#Take the genotypes in any order.
			#REQUIRE_UNIQUE
			#Require that all samples/genotypes be unique between all inputs. 
		]
		code+=map(lambda x: '--variant ' + x, self.infile)
		return '   \\\n\t'.join(code)
	
	
	def horizontalMergeVCFbcftools(self):
		'''
		Merge multiple VCF/BCF files from non-overlapping sample sets to create one multi-sample file. For example, when merging file A.vcf.gz containing samples S1, S2 and S3 and file B.vcf.gz containing samples S3 and S4, the output file will contain four samples named S1, S2, S3, 2:S3 and S4.

		Note that it is responsibility of the user to ensure that the sample names are unique across all files. If they are not, the program will exit with an error unless the option --force-samples is given. The sample names can be also given explicitly using the --print-header and --use-header options.

		Note that only records from different files can be merged, never from the same file. 
		'''
		code=[self.software.bcftools + ' merge',
			'-o ' + self.outfile,  #output FILE
			#'-m '  #The option controls what types of multiallelic records can be created: 
				#-m none   ..  no new multiallelics, output multiple records instead
				#-m snps   ..  allow multiallelic SNP records
				#-m indels ..  allow multiallelic indel records
				#-m both   ..  both SNP and indel records can be multiallelic
				#-m all    ..  SNP records can be merged with indel records
				#-m id     ..  merge by ID
			'-O v'  #--output-type b|u|z|v   Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unecessary compression/decompression and VCF←→BCF conversion. 
			#'-r '  #--regions chr|chr:pos|chr:from-to|chr:from-[,…] 
			#'-R '  #--regions-file file
			#'-l '  #--file-list FILE   read file names from FILE
			#'--force-samples '  #if the merged files contain duplicate samples names, proceed anyway. Duplicate sample names will be resolved by prepending index of the file as it appeared on the command line to the conflicting sample name (see 2:S3 in the above example). 
			#'--print-header'  #print only merged header and exit 
			#' --use-header FILE '  #use the VCF header in the provided text FILE
			#'-f '  #--apply-filters LIST
			#'-i '  #--info-rules -|TAG:METHOD[,…]   Rules for merging INFO fields (scalars or vectors) or - to disable the default rules. METHOD is one of sum, avg, min, max, join. Default is DP:sum,DP4:sum if these fields exist in the input files. Fields with no specified rule will take the value from the first input file. The merged QUAL value is currently set to the maximum. This behaviour is not user controllable at the moment. 
		]
		code+=self.infile
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		if self.type == 'horizontal':
			#subprocess.check_call(self.horizontalMergeVCFbcftools(), shell=True)
			subprocess.check_call(self.horizontalMergeVCFgatk(), shell=True)
		else:
			subprocess.check_call(self.verticalMergeVCFgatk(), shell=True)
			#subprocess.check_call(self.verticalMergeVCFpicard(), shell=True)
		
	def printJob(self):
		if self.type == 'horizontal':
			#return self.horizontalMergeVCFbcftools()
			return self.horizontalMergeVCFgatk() 
		else:
			return self.verticalMergeVCFgatk() 
			#return self.verticalMergeVCFpicard() 
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="vcf files, seperated by ','",required=True)
	parser.add_argument('--outfile', help="merged vcf", required=True)
	parser.add_argument('--type', help="horizontal or vertical merge", default='horizontal', choices=['horizontal','vertical'])
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	mergeVcf=MERGEVCF(args.infile, args.outfile, args.type)
	mergeVcf()