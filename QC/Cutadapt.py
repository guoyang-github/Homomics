# coding=utf-8
# guoyang

import subprocess
import argparse
import os.path
import sys
import os

usage='''
http://cutadapt.readthedocs.io/en/stable/
fqcheck及ng_qc程序，对接头的过滤标准（大致是10个以上碱基的比对，mismatch比率1/10）一些时候并不能完全去除接头，造成GC分布reads末端分叉
需要在这种特殊情况下有一种更完全的接头处理方法

http://genome.annoroad.com/News/Industry/469.html 一个使用文档
'''

class innerSoftware:
	def __init__(self):
		self.java1_8='java'
		self.picard='picard.jar'
		

class BAMINDEX:
	def __init__(self, infile, software=innerSoftware()):
		self.infile=os.path.abspath(infile)
		self.software=software
		
	def BamIndex(self):
		code=[self.software.java1_8 + ' -jar ' + self.software.picard + ' BuildBamIndex',
			'INPUT=' + self.infile
			#'OUTPUT=File'  #The BAM index file. Defaults to x.bai if INPUT is x.bam, otherwise INPUT.bai.
		]		
		return '   \\\n\t'.join(code)
		
		
	def printJob(self):
		return self.BamIndex() + '\n\n'

	def __call__(self):
		subprocess.check_call(self.BamIndex(), shell=True)
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="bam file",required=True)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	bamIndex=BAMINDEX(args.infile)
	bamIndex()