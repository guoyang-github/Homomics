# coding=utf-8
# independent file to encapsulate the software with nesserray format transfor for use independently
import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar picard.jar MarkDuplicates \ 
    INPUT=sorted_reads.bam \ 
    OUTPUT=dedup_reads.bam \
    METRICS_FILE=metrics.txt
'''

class innerSoftware:
	def __init__(self):
		self.java1_8='java'
		self.picard='picard.jar'
		

class MARKDUP:
	'''
	 The tool's main output is a new SAM or BAM file in which duplicates have been identified in the SAM flags field, or optionally removed (see REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES), and optionally marked with a duplicate type in the 'DT' optional attribute. In addition, it also outputs a metrics file containing the numbers of READ_PAIRS_EXAMINED, UNMAPPED_READS, UNPAIRED_READS, UNPAIRED_READ DUPLICATES, READ_PAIR_DUPLICATES, and READ_PAIR_OPTICAL_DUPLICATES. 
	'''
	def __init__(self, infile, outdir, software=innerSoftware()):
		self.infile=infile
		self.outdir=os.path.abspath(outdir)
		if not os.path.isdir(self.outdir):
			os.mkdir(self.outdir)
		self.software=software
		
	def MarkDuplicates(self):
		self.outfile=os.path.join(self.outdir, self.infile.replace('.bam','') + '.dedup.bam')
		self.metrics=self.outfile.replace('.bam','') + '.metrics'
		code=[self.software.java1_8 + ' -Xmx10g -jar ' + self.software.picard + ' MarkDuplicates',
			'INPUT=' + self.infile,   #The BAM or SAM file to sort. Required. 
			'OUTPUT=' + self.outfile,
			'METRICS_FILE=' + self.metrics,  #File to write duplication metrics to Required 
			'TMP_DIR=' + os.path.join(self.outdir,'tmpDeDup')
			#'VALIDATION_STRINGENCY=SILENT',
			#'ASSUME_SORTED=true',
			#'CREATE_INDEX=true',
			#'REMOVE_DUPLICATES=false',
			#'MAX_RECORDS_IN_RAM=4000000'
		]   		
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		subprocess.check_call(self.MarkDuplicates(), shell=True)
		
	def printJob(self):
		return self.MarkDuplicates() + '\n\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="bam file",required=True)
	parser.add_argument('--outdir',help="directory to put output files with suffix 'dedup'", default='.')
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	markDup=MARKDUP(args.infile, args.outdir)
	markDup()	