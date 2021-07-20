# coding=utf-8
# guoyang
# 封装脚本
import subprocess
import argparse
import os.path
import sys
import os

usage='''
	fastqc seqfile1 seqfile2 .. seqfileN -o .

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
           [-c contaminant file] seqfile1 .. seqfileN


Per Tile Sequence Quality

This graph will only appear in your analysis results if you're using an Illumina library which retains its original sequence identifiers. Encoded in these is the flowcell tile from which each read came. The graph allows you to look at the quality scores from each tile across all of your bases to see if there was a loss in quality associated with only one part of the flowcell.
The plot shows the deviation from the average quality for each tile. The colours are on a cold to hot scale, with cold colours being positions where the quality was at or above the average for that base in the run, and hotter colours indicate that a tile had worse qualities than other tiles for that base. In the example below you can see that certain tiles show consistently poor quality. A good plot should be blue all over.



'''

class innerSoftware:
	def __init__(self):
		self.fastqc='fastqc'
		

class FASTQC:
	def __init__(self, infiles, outdir, software=innerSoftware()):
		self.infiles = infile.strip().split(',')
		self.outdir = os.path.absname(outdir)
		self.software=software
		
	def FastQC(self):
		code=[self.software.fastqc,
			'-o ' + self.outdir,
		]
		code += self.infiles
		return '   \\\n\t'.join(code)
		
		
	def printJob(self):
		return self.FastQC() + '\n\n'

	def __call__(self):
		subprocess.check_call(self.printJob(), shell=True)
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infiles',help="fastq files, seperated by ',' ",required=True)
	parser.add_argument('--outdir',help="output dir", default = '.')
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	fastQC=FASTQC(args.infiles, args.outdir)
	fastQC()