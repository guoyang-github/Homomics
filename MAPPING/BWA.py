# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
bwa mem -M -R ’<read group info>’ -p reference.fa raw_reads.fq > aligned_reads.sam
'''

class innerSoftware:
	def __init__(self):
		self.bwa='bwa'
		self.samtools='samtools'

class innerDatabase:
	def __init__(self):
		self.reference='human_g1k_v37_decoy.fasta'
		self.referenceIndex='human_g1k_v37_decoy.fasta.fai'

class BWA:
	def __init__(self, infile, outfile, readgroup=None, software=innerSoftware(), database=innerDatabase()):
		self.inFileList=infile.split(',')
		self.outfile=os.path.abspath(outfile)
		if not os.path.isdir(os.path.dirname(self.outfile)):
			os.mkdir(os.path.dirname(self.outfile))
		self.readgroup=readgroup
		self.software=software
		self.database=database
		
	def bwaMem(self):
		code = [self.software.bwa + ' mem',
			'-M',   #the flag causes BWA to mark shorter split hits as secondary (essential for Picard compatibility)
			'-t 8']
		if self.readgroup:
			code += ['-R ' + self.readgroup] #read group header line such as '@RG\tID:foo\tSM:bar' [null]  \'@RG\\tID:'+self.sampleID+'_'+aLib+'\\tSM:'+self.sampleID+'\\tLB:'+self.sampleID+'\\tPU:'+aLib+'\\tPL:illumina\\tCN:novogene\
		code += self.database.reference
		code += self.inFileList   #fq1,fq2
		code += ['| ' + self.software.samtools + ' view',
			'-b',   #Output in the BAM format
			'-S',   #Ignored for compatibility with previous samtools versions. Previously this option was required if input was in SAM format, but now the correct format is automatically detected by examining the first few characters of input.
			'-t ' + self.database.referenceIndex,   #A tab-delimited FILE. Each line must contain the reference name in the first column and the length of the reference in the second column, with one line for each distinct reference. Any additional fields beyond the second column are ignored. This file also defines the order of the reference sequences in sorting.
			'- > ' + self.outfile]
		
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		subprocess.check_call(self.bwaMem(), shell=True)
		
	def printJob(self):
		return self.bwaMem() + '\n\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="fastq file, seperated by ',' for fq1,fq2",required=True)
	parser.add_argument('--outfile',help="output bam file", required=True)
	parser.add_argument('--readgroup',help="read group header line such as   '@RG\tID:foo\tSM:bar'  '@RG\tID:ZJU2088278T_DHE01682_HF5GGCCXX_L3\tSM:ZJU2088278T\tLB:ZJU2088278T\tPU:DHE01682_HF5GGCCXX_L3\tPL:illumina\tCN:novogene'" ,default=None)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	bwa=BWA(args.infile, args.outfile, args.readgroup)
	bwa()	