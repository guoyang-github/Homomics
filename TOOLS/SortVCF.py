# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar picard.jar SortVcf \
      I=vcf_1.vcf \
      I=vcf_2.vcf \
      O=sorted.vcf
'''

class innerSoftware:
	def __init__(self):
		self.java1_8='java'
		self.picard='picard.jar'
		

class SORTVCF:
	'''
	Sorts one or more VCF files. This tool sorts the records in VCF files according to the order of the contigs in the header/sequence dictionary and then by coordinate. It can accept an external sequence dictionary. If no external dictionary is supplied, the VCF file headers of multiple inputs must have the same sequence dictionaries.

	If running on multiple inputs (originating from e.g. some scatter-gather runs), the input files must contain the same sample names in the same column order. 
	'''
	def __init__(self, infile, outdir, software=innerSoftware()):
		self.infile=os.path.abspath(infile)
		self.outdir=os.path.abspath(outdir)
		self.software=software
		
	def sortVcf(self):
		self.sortedvcf=os.path.join(self.outdir, os.path.basename(self.infile).replace('.vcf','') + '.sorted.vcf')
		code=[self.software.java1_8 + ' -jar ' + self.software.picard + ' SortVcf',
			'INPUT=' + self.infile,
			'OUTPUT=' + self.sortedvcf,
			#'SEQUENCE_DICTIONARY'  #(File)  Default value: null. 
		]		
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		subprocess.check_call(self.sortVcf(), shell=True)
		
	def printJob(self):
		return self.sortVcf() + '\n\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="vcf file",required=True)
	parser.add_argument('--outdir',help="output dir", required=True)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	sortVCF=SORTVCF(args.infile, args.outdir)
	sortVCF()