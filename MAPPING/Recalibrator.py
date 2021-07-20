# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar GenomeAnalysisTK.jar \ 
    -T BaseRecalibrator \ 
    -R reference.fa \ 
    -I realigned_reads.bam \ 
    -L 20 \ 
    -knownSites dbsnp.vcf \ 
    -knownSites gold_indels.vcf \ 
    -o recal_data.table 

java -jar GenomeAnalysisTK.jar \ 
    -T BaseRecalibrator \ 
    -R reference.fa \ 
    -I realigned_reads.bam \ 
    -L 20 \ 
    -knownSites dbsnp.vcf \ 
    -knownSites gold_indels.vcf \ 
    -BQSR recal_data.table \ 
    -o post_recal_data.table 

java -jar GenomeAnalysisTK.jar \ 
    -T AnalyzeCovariates \ 
    -R reference.fa \ 
    -L 20 \ 
    -before recal_data.table \
    -after post_recal_data.table \
    -plots recalibration_plots.pdf

java -jar GenomeAnalysisTK.jar \ 
    -T PrintReads \ 
    -R reference.fa \ 
    -I realigned_reads.bam \ 
    -L 20 \ 
    -BQSR recal_data.table \ 
    -o recal_reads.bam 
'''

class innerSoftware:
	def __init__(self):
		self.java='java'
		self.gatk='GenomeAnalysisTK.jar'
		
class innerDatabase:
	def __init__(self):
		self.reference='human_g1k_v37_decoy.fasta'
		self.Mills1KGindell='Mills_and_1000G_gold_standard.indels.b37.vcf'
		self.dbsnp='00-All.vcf'


class RECALIBRATOR:
	def __init__(self, infile, outdir, software=innerSoftware(), database=innerDatabase()):
		self.infile=infile
		self.prefix=os.path.basename(self.infile).replace('.bam', '')
		self.outdir=os.path.abspath(outdir)
		self.software=software
		self.database=database
		
	def baseRecalibrator(self):
		self.recalTable=os.path.join(self.outdir, self.prefix + '.recal.table'
		code=[self.software.java + ' -Xmx10g -jar ' + self.software.gatk,
			'-T BaseRecalibrator',
			'-nct 6',
			'-R ' + self.database.reference,
			'-I ' + self.infile,
			#'-L 20',
			#'-rf BadCigar',  #filte bad reads
			'-knownSites ' + self.database.dbsnp,
			'-knownSites ' + self.database.Mills1KGindell,
			'-o ' + self.recalTable
		]		
		return '   \\\n\t'.join(code)
		
	def printReads(self):
		self.recalBam=os.path.join(self.outdir, self.prefix + '.recal.bam')
		code=[self.software.java + ' -jar' + self.software.gatk,
			'-T PrintReads',
			'-R' + self.database.reference,
			'-I ' + self.infile,
			'-BQSR ' + self.recalTable,
			'-o ' + self.recalBam
		]
		return '   \\\n\t'.join(code)
		
		
	def __call__(self):
		subprocess.check_call(self.baseRecalibrator(), shell=True)
		subprocess.check_call(self.printReads(), shell=True)
		
	def printJob(self):
		return '\n\n'.join([self.baseRecalibrator(), self.printReads()])
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="bam file",required=True)
	parser.add_argument('--ourdir',help="output dir", default='.')
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	reCal=RECALIBRATOR(args.infile)
	reCal()