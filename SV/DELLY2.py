# coding=utf-8
# independent file to encapsulate the software with nesserray format transfor for use independently
import subprocess
import argparse
import os.path
import sys
import os

usage='''
Delly2 needs a sorted and indexed bam file for every sample and the reference genome to identify split-reads. The output is in bcf format with a csi index. Prior duplicate marking is recommended. The SV type can be DEL, DUP, INV, TRA, or INS for deletions, tandem duplications, inversions, translocations and small insertions, respectively. Delly2 supports germline and somatic SV discovery, genotyping and filtering. Because of that, Delly2 has been modularized and common workflows for germline and somatic SV calling are outlined below. If you do need VCF output you need a recent version of bcftools (included as a submodule in Delly2) for file conversion.

SV calling is done by sample or in small batches to increase SV sensitivity & breakpoint precision.

somatic  germline
DEL, DUP, INV, TRA, or INS
deletions, tandem duplications, inversions, translocations and small insertions


和群体call类似 (breakdancer,GATK) 互做参照   这种方法应该可以做 denovoSV 新生突变分析
Somatic SV calling

At least one tumor sample and a matched control sample are required. All tumor/control pairs are run separately for SV discovery:

	delly call -t DEL -x hg19.excl -o t1.bcf -g hg19.fa tumor1.bam control1.bam

Somatic pre-filtering of every tumor/control pair using a tab-delimited sample description file where the first column is the sample id (as in the VCF/BCF file) and the second column is either tumor or control.

	delly filter -t DEL -f somatic -o t1.pre.bcf -s samples.tsv -g hg19.fa t1.bcf

Re-genotype somatic sites across a larger panel of control samples to efficiently filter false postives and germline SVs. For performance reasons, this can be run in parallel for each sample (see germline SV calling) and/or directly on a combined pre-filtered somatic site list from multiple tumor/control pairs.

	delly call -t DEL -g hg19.fa -v t1.pre.bcf -o geno.bcf -x hg19.excl tumor1.bam control1.bam ... controlN.bam   是不是并不改变－v文件中的位置，只是找了下后面文件的基因型，可能只是过滤假阳性吧    有位点参照（-v）的重新call一下，过滤可能的假阳性somatic

Post-filter for somatic SVs using all control samples.

	delly filter -t DEL -f somatic -o t1.somatic.bcf -s samples.tsv -g hg19.fa geno.bcf
	
	

Germline SV calling

SV calling is done by sample or in small batches to increase SV sensitivity & breakpoint precision

	delly call -t DEL -g hg19.fa -o s1.bcf -x hg19.excl sample1.bam

Merge SV sites into a unified site list

	delly merge -t DEL -m 500 -n 1000000 -o del.bcf -b 500 -r 0.5 s1.bcf s2.bcf ... sN.bcf     联合校正call       这个后面的merged.bcf有啥区别？

Re-genotype merged SV site list across all samples. This can be run in parallel for each sample.

	delly call -t DEL -g hg19.fa -v del.bcf -o s1.geno.bcf -x hg19.excl s1.bam          

	delly call -t DEL -g hg19.fa -v del.bcf -o sN.geno.bcf -x hg19.excl sN.bam

Merge all re-genotyped samples to get a single VCF/BCF using bcftools merge

	bcftools merge -O b -o merged.bcf s1.geno.bcf s2.geno.bcf ... sN.geno.bcf

Apply the germline SV filter

	delly filter -t DEL -f germline -o germline.bcf -g hg19.fa merged.bcf



What is the smallest SV size Delly can call?
This depends on the sharpness of the insert size distribution. For an insert size of 200-300bp with a 20-30bp standard deviation, Delly starts to call reliable SVs >=300bp. Delly2 also supports calling of small InDels using soft-clipped reads only. In this mode the smallest SV size called is 15bp.



1. For each input BAM file, DELLY computes the default read-pair orientation and the paired-end insert size distribution characterized by the median and standard deviation of the library. Based on these parameters, DELLY then identifies all discordantly mapped read-pairs that either have an abnormal orientation or an insert size greater than the expected range. 
2. The paired-end clusters identified in the previous mapping analysis are interpreted as breakpoint-containing genomic intervals, which are subsequently screened for split-read support to fine map the genomic rearrangements at single-nucleotide resolution and to investigate the breakpoints for potential microhomologies and microinsertions.

DELLY has been specifically geared towards enabling SV calling in the presence of different paired-end sequencing libraries with distinct insert sizes

做个体样本的SV calling，有关群体矫正的步骤全部不做，不适合流程化
'''

class innerSoftware:
	def __init__(self):
		self.delly2='delly_v0.7.3_parallel_linux_x86_64bit'
		self.bcftools='bcftools'

class innerDatabase:
	def __init__(self):
		self.reference = 'human_g1k_v37_decoy.fasta'
		#self.excl= ''
		

class DELLY:
	def __init__(self, case, control, outdir, type, software=innerSoftware(), database=innerDatabase()):
		self.case=case
		self.control=control
		self.outdir=os.path.abspath(outdir)
		self.type=type
		self.software=software
		self.database=database
	
	def sampleDesc(self):
		self.sampleDesc=os.path.join(self.outdir, 'sampledescription.txt')
		desc=[]
		desc.append(os.path.basename(self.case).replace('.bam','') + '\ttumor\n')
		desc.append(os.path.basename(self.control).replace('.bam','') + '\tcontrol\n')
		with open(self.sampleDesc,'w') as sampleDescFile:
			sampleDescFile.writelines(desc)
	
		
	def germline(self):
		callCode=[self.software.delly2 + ' call ',
			'-t ' + self.type,
			'-g ' + self.database.reference,
			'-q 20',  #min. paired-end mapping quality
			'-o ' + os.path.join(self.outdir, os.path.basename(self.case) + '.bcf'),
			#'-x ' + self.database.excl,   #Exclude telomere and centromere regions
			'--noindels', #no small InDel calling
			#'--indelsize 500', #max. small InDel size
			self.case
		]
		self.germlineBcf=os.path.join(self.outdir, os.path.basename(self.case) + '.germline.bcf')
		filterCode=[self.software.delly2 + ' filter ', 
			'-t ' + self.type,
			'-f germline',
			'-o ' + self.germlineBcf,
			'-g ' + self.database.reference,
			os.path.join(self.outdir, os.path.basename(self.case) + '.bcf')
		]
		self.germlineVcf=os.path.join(self.outdir, os.path.basename(self.case) + '.germline.vcf')
		transCode=self.software.bcftools + ' view ' + self.germlineBcf + ' > ' + self.germlineVcf	
		return '   \\\n\t'.join(callCode) + '\n' + '   \\\n\t'.join(filterCode) + '\n' + transCode
		
	def somatic(self):
		#从somatic mutations 中用群体样本集to efficiently filter false postives and germline SVs 不做并且做也不一定合适
		callCode=[self.software.delly2 + ' call ',
			'-t ' + self.type,
			#'-x ' + self.database.excl,
			'-o ' + os.path.join(self.outdir, os.path.basename(self.case) + '.bcf'),
			'-g ' + self.database.reference,
			'-q 20',  #min. paired-end mapping quality
			'--noindels', #no small InDel calling
			#'--indelsize 500', #max. small InDel size
			self.case,
			self.control			
		]
		self.somaticBcf=os.path.join(self.outdir, os.path.basename(self.case) + '.somatic.bcf')
		filterCode=[self.software.delly2+ ' filter '，
			'-t ' + self.type,
			'-f somatic',
			'-o ' + self.germlineBcf,
			'-g ' + self.database.reference,
			'-s ' + self.sampleDesc,
			os.path.join(self.outdir, os.path.basename(self.case) + '.bcf')
		]
		self.somaticVcf=os.path.join(self.outdir, os.path.basename(self.case) + '.somatic.vcf')
		transCode=self.software.bcftools + ' view ' + self.somaticBcf + ' > ' + self.somaticVcf
		return '   \\\n\t'.join(callCode) + '\n' + '   \\\n\t'.join(filterCode) + '\n' + transCode
		
	def __call__(self):
		if self.control:
			self.sampleDesc()
			subprocess.check_call(self.somatic(), shell=True)
		else:
			subprocess.check_call(self.germline(), shell=True)
		
	def printJob(self):
		if self.control:
			self.sampleDesc()
			return self.somatic() + '\n\n'
		else:
			return self.germline() + '\n\n'


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--case', help="case bam file",required=True)
	parser.add_argument('--control', help="control bam file for somatic", default=None)
	parser.add_argument('--outdir', help="output dir", default='.')
	parser.add_argument('--type', help="SV type: DEL, DUP, INV, TRA, or INS", choices=['DEL', 'DUP', 'INV', 'TRA', 'INS'], required=True)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	delly=DELLY(args.case, args.control, args.outdir, args.type)
	delly()