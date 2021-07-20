# coding=utf-8
# independent file to encapsulate the software with nesserray format transfor for use independently
import subprocess
import argparse
import os.path
import sys
import shutil
import os

usage='''
ANNOVAR input format
The first five space or tab delimited fields are Chromosome ("chr" prefix is optional), Start, End, Reference Allelel, Alternative Allele. 
The rest of the columns are completely optional.


usage: python Annovar.py -in example.vcf -build hg19 -out out/example
-in: vcf file or annovar input format, example: /path/to/input/example.vcf or /path/to/input/example.avinput
-out: directory to put output files, example: path/to/out/
-build: build version, default: hg19


table_annovar.pl 可以接受vcf，并输出vcf ／txt             流程化的，输出所有类型的注释信息，变成一张大表 (一般都用这个脚本)
	In addition to ANNOVAR input files, table_annovar can directly take VCF file as input now, 
	and generate a VCF file (as well as tab-delimited text file) as output file, 
	with its INFO field populated with various ANNOVAR annotations.
annotate_variation.pl 需要annovar格式的输入              分步骤的，annovar基础程序  下载数据  分 gene region filter
convert2annovar.pl  可以把其他格式的文件（vcf等）   变成annovar格式


So as ANNOVAR developer, I decided to re-process all 1000 Genomes Project files as well as ESP6500si files as well as dbSNP files, 
so that each line contains one variant and so that every variant is left-normalized. 
The updated databases were made available in December 2014.
So as a user, this is what you should do: 
(1) split VCF lines so that each line contains one and only one variant 
(2) left-normalize all VCF lines 
(3) annotate by ANNOVAR.
HGVS clearly specifies that left-normalization would be performed on cDNA (mRNA) coordinate, which means that right-normalization is required for half the genes in human genome.



Gene-based annotation: identify whether SNPs or CNVs cause protein coding changes and the amino acids that are affected. Users can flexibly use RefSeq genes, UCSC genes, ENSEMBL genes, GENCODE genes, AceView genes, or many other gene definition systems.
Region-based annotation: identify variants in specific genomic regions, for example, conserved regions among 44 species, predicted transcription factor binding sites, segmental duplication regions, GWAS hits, database of genomic variants, DNAse I hypersensitivity sites, ENCODE H3K4Me1/H3K4Me3/H3K27Ac/CTCF sites, ChIP-Seq peaks, RNA-Seq peaks, or many other annotations on genomic intervals.
Filter-based annotation: identify variants that are documented in specific databases, for example, whether a variant is reported in dbSNP, what is the allele frequency in the 1000 Genome Project, NHLBI-ESP 6500 exomes or Exome Aggregation Consortium, calculate the SIFT/PolyPhen/LRT/MutationTaster/MutationAssessor/FATHMM/MetaSVM/MetaLR scores, find intergenic variants with GERP++ score < 2, or many other annotations on specific mutations.


'''

class innerSoftware:
	def __init__(self):
		self.tableAnnovar='table_annovar.v2.pl'
		self.bcftools='bcftools'
		self.bgzip='bgzip'
		self.tabix='tabix'
		
class innerDatabase:
	def __init__(self):
		self.anndatabase='humandb_b37/'
		self.reference='human_g1k_v37_decoy.fasta'

class ANNOVAR:
	def __init__(self, infile, outdir, build='hg19', software=innerSoftware(), database=innerDatabase()):
		self.infile=os.path.abspath(infile)
		self.infileBasename=os.path.basename(self.infile)
		self.outdir=os.path.abspath(outdir)
		if not os.path.isdir(self.outdir):
			os.mkdir(self.outdir)
		self.build=build
		self.software=software
		self.database=database
	
	def gbZipIndex(self):
		if not os.path.samefile(os.path.dirname(self.infile),self.outdir):
			shutil.copy(self.infile, self.outdir)
		code = []
		if not self.infile.endswith('.gz'):
			codeBgzip=self.software.bgzip + ' ' + os.path.join(self.outdir,self.infileBasename)
			code.append(codeBgzip)
		self.zipInFile=os.path.join(self.outdir,self.infileBasename + '.gz')
		self.zipInFileIndex = self.zipInFile + '.tbi'
		if not os.path.exists(self.zipInFileIndex):
			codeTabix=self.software.tabix + ' ' + self.zipInFile
			code.append(codeTabix)
		return '\n'.join(code)
						
	
	def splitLine(self):
		self.splitLineInput=os.path.join(self.outdir, self.infileBasename.replace('.vcf','') + '.splitLine.vcf') 
		code=[self.software.bcftools + ' norm',
			'-m-both ',
			'-o ' + self.splitLineInput,
			self.zipInFile]
		return '   \\\n\t'.join(code)
		
	def lifeNorm(self):
		self.splitLineLifeNormInput=self.splitLineInput.replace('.vcf','')+'.lifeNorm.vcf'
		code=[self.software.bcftools + ' norm',
			'-f ' + self.database.reference,
			'-o ' + self.splitLineLifeNormInput,
			self.splitLineInput]
		return '   \\\n\t'.join(code)
		
	def annoVar(self):
		code=[self.software.tableAnnovar,
			self.splitLineLifeNormInput,
			self.database.anndatabase,
			'-buildver ' + self.build,
			'-out ' + os.path.join(self.outdir, self.infileBasename.replace('.vcf','')),  #output file name prefix
			'-remove ',
			'-protocol refGene,1000g2014aug_all,snp138,cosmic70,clinvar_20150629,esp6500si_all,ljb26_all,SAOdbSNP141,SAOclinvar20140303,exac03',
			'-operation g,f,f,f,f,f,f,f,f,f ',
			"--argument '--splicing_threshold 10 --transcript_function,,,,,,,,,'",
			'-nastring . ',
			'-vcfinput']
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		subprocess.check_call(self.printJob(), shell=True)

	
	def printJob(self):
		return '\n\n'.join([self.gbZipIndex(), self.splitLine(), self.lifeNorm(), self.annoVar()])+'\n'
		

def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="vcf file",required=True)
	parser.add_argument('--outdir',help="directory to put output files",default='.')
	parser.add_argument('--build',help="build version",default='hg19', choices=['hg19','b38'])
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	annoVar=ANNOVAR(args.infile, args.outdir, args.build)
	annoVar()

