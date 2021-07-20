# coding=utf-8
# guoyang

import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar picard.jar BuildBamIndex \ 
    INPUT=dedup_reads.bam
'''

class innerSoftware:
	def __init__(self):
		self.java1_8 = 'java'
		self.picard = 'picard.jar'
		

class BAMINDEX:
	def __init__(self, infile, software = innerSoftware()):
		self.infile = os.path.abspath(infile)
		self.software = software
		
	def BamIndex(self):
		code = [self.software.java1_8 + ' -jar ' + self.software.picard + ' BuildBamIndex',
			'INPUT=' + self.infile
			#'OUTPUT=File'  #The BAM index file. Defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai.
		]		
		return '   \\\n\t'.join(code)
		
		
	def printJob(self):
		return self.BamIndex() + '\n\n'

	def __call__(self):
		subprocess.check_call(self.BamIndex(), shell = True)
		


def paramsParse():
	parser = argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument('--infile', help = "bam file", required = True)
	return parser.parse_args()

if __name__ == '__main__':
	args = paramsParse()
	bamIndex = BAMINDEX(args.infile)
	bamIndex()