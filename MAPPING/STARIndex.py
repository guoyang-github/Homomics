# coding=utf-8
# guoyang
import subprocess
import argparse
import os.path
import sys
import os

usage='''
STAR
--runThreadN NumberOfThreads
--runMode genomeGenerate
--genomeDir /path/to/genomeDir
--genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ... --sjdbGTFfile /path/to/annotations.gtf
--sjdbOverhang ReadLength-1



参考序列的选择
Which chromosomes/scaffolds/patches to include?
It is strongly recommended to include major chromosomes (e.g., for human chr1-22,chrX,chrY,chrM,) 
as well as un-placed and un-localized scaffolds. Typically, un-placed/un-localized scaffolds add just 
a few MegaBases to the genome length, however, a substantial number of reads may map to ribosomal RNA (rRNA) repeats 
on these scaffolds. These reads would be reported as unmapped if the scaffolds are not included in the genome, 
or, even worse, may be aligned to wrong loci on the chromosomes. 
Generally, patches and alternative haplotypes should not be included in the genome.
'''

class innerSoftware:
	def __init__(self):
		self.star = ''

		

class STARINDEX:
	def __init__(self, runThreadN, genomeDir, genomeFastaFiles, sjdbGTFfile, sjdbOverhang, annotFormat, software = innerSoftware()):
		self.runThreadN = runThreadN
		self.genomeDir = genomeDir   #结果文件夹，可当成对象属性引出给下游对象使用
		if not os.path.exist(self.genomeDir):
			os.mkdir(self.genomeDir)
		self.genomeFastaFiles = genomeFastaFiles
		self.sjdbGTFfile = sjdbGTFfile
		self.sjdbOverhang = sjdbOverhang
		self.annotFormat = annotFormat
		self.software = software
		
	def StarIndex(self):
		code = [self.software.star,
			'--runThreadN ' + self.runThreadN,
			'--genomeDir ' + self.genomeDir,
			'--genomeFastaFiles' + self.genomeFastaFiles,
			'--sjdbGTFfile ' + self.sjdbGTFfile,
			'--sjdbOverhang ' + self.sjdbOverhang
		]	
		if self.annotFormat == 'GFF':
			code += ['--sjdbGTFtagExonParentTranscript Parent']	
		return '   \\\n\t'.join(code)
		
		
	def printJob(self):
		return self.StarIndex() + '\n\n'

	def __call__(self):
		subprocess.check_call(self.StarIndex(), shell = True)
		


def paramsParse():
	parser = argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)    
	parser.add_argument('--runThreadN', help = "NumberOfThreads", default = 8)
	parser.add_argument('--genomeDir', help = "specifies path to the genome directory", required = True)
	parser.add_argument('--genomeFastaFiles', help = "specifies one or more FASTA files with the genome reference sequences", required = True)
	parser.add_argument('--sjdbGTFfile', help = "specifies the path to the file with annotated transcripts in the standard GTF format")
	parser.add_argument('--sjdbOverhang', help = "ReadLength-1", default = 150 - 1)
	parser.add_argument('--annotFormat', help = 'GTF/GFF3', default = 'GTF', choices = ['GTF', 'GFF'])
	return parser.parse_args()

if __name__ == '__main__':
	args = paramsParse()
	starIndex = STARINDEX(args.runThreadN, args.genomeDir, args.genomeFastaFiles, args.sjdbGTFfile, args.sjdbOverhang, args.annotFormat)
	starIndex()