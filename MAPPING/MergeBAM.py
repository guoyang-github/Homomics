# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar picard.jar MergeSamFiles \
      I=input_1.bam \
      I=input_2.bam \
      O=merged_files.bam
'''

class innerSoftware:
	def __init__(self):
		self.java1_8='java'
		self.picard='picard.jar'
		

class MERGEBAM:
	def __init__(self, infile, outfile, software=innerSoftware()):
		self.infile=infile.strip().split(',')
		self.outfile=os.path.abspath(outfile)
		self.software=software
		
	def mergeBAM(self):
		assert len(self.infile) != 1
		code=[self.software.java1_8 + ' -jar ' + self.software.picard + ' MergeSamFiles']
		code+=map(lambda x: 'INPUT=' + x, self.infile)
		code+=[
			'OUTPUT=' + self.outfile,
			'SORT_ORDER=coordinate',
			#'TMP_DIR=' + os.path.join(oa.path.dirname(self.outfile),'tmp')
			'USE_THREADING=true',
			#'CREATE_INDEX=true',
			#'VALIDATION_STRINGENCY=SILENT',
			#'MAX_RECORDS_IN_RAM=4000000',
			#'MERGE_SEQUENCE_DICTIONARIES',  #Merge the sequence dictionaries Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false} 	
		]
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		subprocess.check_call(self.mergeBAM(), shell=True)
		
	def printJob(self):
		return self.mergeBAM() + '\n\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="bam file list, seperated by ','",required=True)
	parser.add_argument('--outfile', help="merged bam", required=True)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	mergeBam=MERGEBAM(args.infile)
	mergeBam()