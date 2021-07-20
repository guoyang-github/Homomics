# coding=utf-8
# guoyang
import subprocess
import argparse
import os.path
import sys
import os

usage='''
Spliced Transcripts Alignment to a Reference

'''

class innerSoftware:
	def __init__(self):
		self.star = ''
		
		

class STAR:
	def __init__(self, runThreadN, genomeDir, readFile, outFileNamePrefix, software = innerSoftware()):
		self.runThreadN = runThreadN
		self.genomeDir = genomeDir
		self.readFile = ' '.join(readFile.split(','))
		self.outFileNamePrefix = outFileNamePrefix
		
	def Star(self):
		code = [self.software.star,
			'--runThreadN ' + self.runThreadN,
			'--genomeDir ' + self.genomeDir,
			'--readFilesIn ' + self.readFile,  #paths to files that contain input read1 (and, if needed, read2)
			'--outFileNamePrefix ' + self.outFileNamePrefix,  #output files name prefix (including full or relative path)
			'--outSAMtype BAM',
			'--outSAMstrandField intronMotif',   #For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute, which STAR will generate with --outSAMstrandField intronMotif option. As required, the XS strand attribute will be generated for all alignments that contain splice junctions. The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed.
			'--outSAMattrIHstart 0'   #default: 1  start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie
		]		
		return '   \\\n\t'.join(code)
		
		
	def printJob(self):
		return self.Star() + '\n\n'

	def __call__(self):
		subprocess.check_call(self.printJob(), shell = True)
		


def paramsParse():
	parser = argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument('--runThreadN', help = "NumberOfThreads", default = 8)
	parser.add_argument('--genomeDir', help = "specifies path to the genome directory", required = True)
	parser.add_argument('--readFile', help = "read file(s), seperated by ','", required = True)
	parser.add_argument('--outFileNamePrefix', help = "output file prefix", default = './')

	return parser.parse_args()

if __name__ == '__main__':
	args = paramsParse()
	star = STAR(args.runThreadN, args.genomeDir, args.readFile, args.outFileNamePrefix)
	star()