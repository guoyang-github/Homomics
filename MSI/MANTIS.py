# coding=utf-8
# guoyang

import subprocess
import argparse
import os.path
import sys
import os

usage='''
1. 通过RepeatFinder先找MS loci （for finding microsatellites within a reference genome） 看来是什么物种都可以了
2. 已知loci集之后，对每个位置计算normal和tumor的差别stepwise difference 公式很简单，见文献
3. 把所有位点的值算个平均
4. 通过已知数据集计算出一个最好的阈值0.4，并以此为依据

根据不同测序策略的测序深度，在filter的过程中阈值会稍有不同，软件默认是按照target测序来的，本程序中改成按照WES的来
'''

class innerSoftware:
	def __init__(self):
		self.repeatfinder = ''
		self.mantis = ''
		self.python = ''
		

class MANTIS:
	def __init__(self, normal, tumor, bedfile, genome, output, software = innerSoftware()):
		self.normal = os.path.abspath(normal)
		self.tumor = os.path.abspath(tumor)
		self.bedfile = bedfile
		self.genome = os.path.abspath(genome)
		self.output = os.path.abspath(output)
		self.software = software
		
	def RepeatFinder(self):
		self.bedfile = os.path.join(os.path.dirname(self.output), self.genome + '.ms.loci.bed')
		code = [self.software.repeatfinder,
			'-i ' + self.genome,
			'-o ' + self.bedfile
		]		
		return '   \\\n\t'.join(code)

	def Mantis(self):
		code = [self.software.python + ' ' + self.software.mantis,
			'--bedfile ' + self.bedfile,
			'--genome ' + self.genome,
			'-n ' + self.normal,
			'-t ' + self.tumor,
			'-o ' + self.output,
			#default quality thresholds are intended for use in situations such as targeted resequencing, 
			#in which locus coverage is less of an issue than with whole-exome data.
			#a less stringent set of thresholds for whole-exome data
			'-mrq 20.0', 
			'-mlq 25.0', 
			'-mlc 20',
			'-mrr 1'
		]
		return '   \\\n\t'.join(code)
		
	def printJob(self):
		if self.bedfile:
			return self.Mantis() + '\n\n'
		else:
			return '\n\n'.join(self.RepeatFinder(), self.Mantis())

	def __call__(self):
		subprocess.check_call(self.printJob(), shell = True)
		


def paramsParse():
	parser = argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument('--normal', help = "normal bam file", required = True)
	parser.add_argument('--tumor', help = "tumor bam file", required = True)
	parser.add_argument('--bedfile', help = "customized bed file", default = None)
	parser.add_argument('--genome', help = "reference genome", required = True)
	parser.add_argument('--output', help = "output file", required = True)

	return parser.parse_args()

if __name__ == '__main__':
	args = paramsParse()
	mantis = MANTIS(args.normal, args.tumor, args.bedfile, args.genome, args.output)
	mantis()