# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar picard.jar SortSam \ 
    INPUT=aligned_reads.sam \ 
    OUTPUT=sorted_reads.bam \ 
    SORT_ORDER=coordinate
'''

class innerSoftware:
	def __init__(self):
		self.java1_8='java'
		self.picard='picard.jar'
		

class SORTBAM:
	def __init__(self, infile, outdir, order, software=innerSoftware()):
		self.infile=infile
		self.outdir=os.path.abspath(outdir)
		self.outfile=os.path.join(self.outdir, self.infile.replace('.bam','') + '.sorted.bam')
		if not os.path.isdir(self.outdir):
			os.mkdir(self.outdir)
		self.order = order
		self.software=software
		
	def SortSam(self):
		code=[self.software.java1_8 + ' -jar ' + self.software.picard + ' SortSam',
			'INPUT=' + self.infile,   #The BAM or SAM file to sort. Required. 
			'OUTPUT=' + self.outfile,
			'SORT_ORDER=' + self.order]   #Sort order of output file Required. Possible values: {unsorted, queryname, coordinate, duplicate} 		
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		subprocess.check_call(self.SortSam(), shell=True)
		
	def printJob(self):
		return self.SortSam() + '\n\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="bam file",required=True)
	parser.add_argument('--outdir',help="directory to put output files with suffix 'sorted'", default='.')
	parser.add_argument('--order',help="Sort order [unsorted, queryname, coordinate, duplicate] default: coordinate", default='coordinate', choices = ['unsorted', 'queryname', 'coordinate', 'duplicate'])
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	sortBam=SORTBAM(args.infile, args.outdir, args.order)
	sortBam()	