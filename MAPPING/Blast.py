# coding=utf-8
# guoyang
# 未完成
import subprocess
import argparse
import os.path
import sys
import os

usage='''
Download pre-formatted BLAST databases from NCBI
update_blastdb.pl [options] blastdb ...

--showall
	Show all available pre-formatted BLAST databases (default: false). The output of this option lists the database names
    which should be used when requesting downloads or updates using this script.



格式化数据库
makeblastdb -in db.fasta -dbtype prot -parse_seqids -out dbname
参数说明:
-in：待格式化的序列文件
-dbtype：数据库类型，prot或nucl
-out：数据库名


蛋白序列比对蛋白数据库（blastp）
blastp -query seq.fasta -out seq.blast -db dbname -outfmt 6 -evalue 1e-5 -num_descriptions 10 -num_threads 8
参数说明:
-query： 输入文件路径及文件名
-out：输出文件路径及文件名
-db：格式化了的数据库路径及数据库名
-outfmt：输出文件格式，总共有12种格式，6是tabular格式对应BLAST的m8格式
-evalue：设置输出结果的e-value值
-num_descriptions：tabular格式输出结果的条数
-num_threads：线程数


核酸序列比对核酸数据库（blastn）以及核酸序列比对蛋白数据库（blastx）
与上面的blastp用法类似：
blastn -query seq.fasta -out seq.blast -db dbname -outfmt 6 -evalue 1e-5 -num_descriptions 10 -num_threads 8
blastx -query seq.fasta -out seq.blast -db dbname -outfmt 6 -evalue 1e-5 -num_descriptions 10 -num_threads 8

以上的参数说明只是一些常用的参数，完整的参数说明可以用-help查询。

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