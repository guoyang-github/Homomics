# coding=utf-8
# 已完成
import subprocess
import argparse
import os.path
import sys
import os

usage='''

It focuses primarily on variants in copy-number neutral, loss of heterozygosity (LOH)-free portions of the genome, 
which allows for the highest-confidence quantification of variant allele frequencies (VAF) and inference of clonality.


library(sciClone)

#read in vaf data from three related tumors
#format is 5 column, tab delimited: 
#chr, pos, ref_reads, var_reads, vaf

v1 = read.table("data/vafs.tumor1.dat",header=T);
v2 = read.table("data/vafs.tumor2.dat",header=T);
v3 = read.table("data/vafs.tumor3.dat",header=T);

#read in regions to exclude (commonly LOH)
#format is 3-col bed
regions = read.table("data/exclude.loh")

#read in segmented copy number data
#4 columns - chr, start, stop, segment_mean   
cn1 = read.table("data/copy_number_tum1")
cn2 = read.table("data/copy_number_tum2")
cn3 = read.table("data/copy_number_tum3")

#set sample names
names = c("Sample1","Sample2","Sample3")


#Examples:
#------------------------------------
#1d clustering on just one sample
sc = sciClone(vafs=v1,
         copyNumberCalls=cn1,
         sampleNames=names[1],
         regionsToExclude=reg1)
#create output
writeClusterTable(sc, "results/clusters1")
sc.plot1d(sc,"results/clusters1.1d.pdf")

#------------------------------------
#2d clustering using two samples:
sc = sciClone(vafs=list(v1,v2),
              copyNumberCalls=list(cn1,cn2),
              sampleNames=names[1:2],
               regionsToExclude=regions)
#create output
writeClusterTable(sc, "results/clusters2")
sc.plot1d(sc,"results/clusters2.1d.pdf")
sc.plot2d(sc,"results/clusters2.2d.pdf")


#------------------------------------
#3d clustering using three samples:
sc = sciClone(vafs=list(v1,v2,v3),
              copyNumberCalls=list(cn1,cn2,cn3),
              sampleNames=names[1:3],
               regionsToExclude=regions)
#create output
writeClusterTable(sc, "results/clusters2")
sc.plot1d(sc,"results/clusters2.1d.pdf")
sc.plot2d(sc,"results/clusters2.2d.pdf")
sc.plot3d(sc, sc@sampleNames, size=700, outputFile="results/clusters3.3d.gif")

#This pattern generalizes up to N samples, except for plotting, which caps out at 3d for obvious reasons.

https://github.com/genome/sciclone


主要的工作就是在按照somatic突变的点的频率来分类，而且只考虑CN2区域的点
另外如果增加时序的取样数据及组织多点取样的样本的话，可以统筹考虑，把分类做的更准更细
比如两种突变的频率差不多，如果只有一个样本，那就会分在一起，如果多一个样本做参考，能看出两组突变的变化差异，就能更好的分组
'''

class innerSoftware:
	def __init__(self):
		self.R3_2 = 'Rscript'
		self.bedtools = 'bedtools'


		

class SCICLONE:
	def __init__(self, vaf, exclude, cn, sample, outdir, software=innerSoftware()):
		self.vaf = vaf
		volume = len(self.vaf.strip().split(','))
		self.exclude = exclude  #exclude有几个都行，最后把这些LOH的区域bed merge成一个
		self.cn = cn
		assert len(self.cn.strip().split(',')) == volume
		self.sample = sample
		assert len(self.sample.strip().split(',')) == volume
		self.outdir = os.path.abspath(outdir)
		self.software = software
		
	
	def mergeBED(self):
		if not self.exclude:
			return ''
		self.excludeCat = os.path.join(self.outdir, 'exclude.regions')
		codeCat = 'cat ' + self.exclude.replace(',',' ') + '> ' self.excludeCat
		self.excludeMerge = os.path.join(self.outdir, 'excludeMerge.regions')
		codeMerge = self.software.bedtools + ' merge -i ' + self.excludeCat + ' > ' + self.excludeMerge
		return codeCat + '\n' + codeMerge + '\n'
			
		
	def sciClone(self):
		self.software.sciClone = os.path.join(sys.path[0],'sciClone.R')
		code = [self.software.R3_2, self.software.sciClone, self.vaf, self.excludeMerge, self.cn, self.sample, self.outdir]	
		return ' '.join(code) + '\n'
		
	def __call__(self):
		return self.mergeBED() + self.sciClone()
		
	def printJob(self):
		subprocess.check_call(self.mergeBED(), shell=True)
		subprocess.check_call(self.sciClone(), shell=True)
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--vaf',help="vaf data, formatted in 5 column(chr, pos, ref_reads, var_reads, vaf), seperated by ','", required=True)
	parser.add_argument('--exclude', help="regions to exclude (commonly LOH), format is 3-col bed", default=None)
	parser.add_argument('--cn', help="read in segmented copy number data, 4 columns - chr, start, stop, segment_mean", required=True)
	parser.add_argument('--sample', help="sample names", required=True)
	parser.add_argument('--outdir', help="output dir", default='.')
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	sciClone = SCICLONE(args.vaf, args.exclude, args.cn, args.samples, args.outdir)
	sciClone()