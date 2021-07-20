# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar GenomeAnalysisTK.jar \
   -T FlagStat \
   -R reference.fasta \
   -I reads.bam \
   [-o output.txt]
'''

class innerSoftware:
	def __init__(self):
		self.gatk='GenomeAnalysisTK.jar'
		self.python='python'
		self.java1_8='java'

class innerDatabase:
	def __init__(self):
		self.reference='human_g1k_v37_decoy.fasta'		

class ALIGNSTAT:
	def __init__(self, infile, outdir, software=innerSoftware(), database=innerDatabase()):
		self.infile=os.path.abspath(infile)
		self.outdir=os.path.abspath(outdir)
		self.software=software
		self.database=database
		
	def AlnStat(self):
		self.software.alnstat = os.path.join(sys.path[0], 'AlnStat.py')
		self.flagstat=os.path.join(self.outdir, os.path.basename(self.infile).replace('.bam','') + '.flagstat')
		code=[self.software.python + ' ' + self.software.alnstat,
			'--bam ' + self.infile,
			'--outdir ' + self.outdir
		]		
		return '   \\\n\t'.join(code)
		
	def AlnStatgatk(self):
		self.flagstatgatk=os.path.join(self.outdir, os.path.basename(self.infile).replace('.bam','') + '.flagstat')
		code=[self.software.java1_8 + ' -jar ' + self.software.gatk,
			'-T FlagStat',
			'-R ' + self.database.reference,
			'-I ' + self.infile,
			'-o ' + self.flagstatgatk
		]
		return '   \\\n\t'.join(code)
		
		
	def __call__(self):
		subprocess.check_call(self.AlnStat(), shell=True)
		
	def printJob(self):
		return self.AlnStat() + '\n\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="bam file",required=True)
	parser.add_argument('--outdir',help="output dir", required=True)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	alnStat=ALIGNSTAT(args.infile, args.outdir)
	alnStat()