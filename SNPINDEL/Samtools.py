# coding=utf-8
# independent file to encapsulate the software with nesserray format transfor for use independently
import subprocess
import argparse
import os.path
import sys
import shutil
import os

usage='''
samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf  
bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf

Variant calling is basically a three-step process:

First, samtools mpileup command transposes the mapped data in a sorted BAM file fully to genome-centric coordinates. 
It starts at the first base on the first chromosome for which there is coverage and prints out one line per base. 
Each line has information on every base observed in the raw data at that base position 
along with a lot of auxiliary information depending on which flags are set. 
It calculates the Bayseian prior probability given by the data, but does not estimate a real genotype.
Next, bcftools with a few options added uses the prior probability distribution and the data to calculate an actual genotype 
for the variants detected.
Finally, vcfutils.pl (or equivalent) is used to filter down the list of candidates according to some set of objective criteria.
'''

class innerSoftware:
	def __init__(self):
		self.reference='human_g1k_v37_decoy.fasta'
		self.samtools='samtools'
		self.bcftools='bcftools'

class SAMTOOLS:
	def __init__(self, bam, outdir, target=None, chrom=None, rawvcf=False, mpileup=False, software=innerSoftware()):
		self.bam=os.path.abspath(bam)
		self.bamBasename=os.path.basename(self.bam)
		self.target=target
		self.outdir=os.path.abspath(outdir)
		if not os.path.isdir(self.outdir):
			os.mkdir(self.outdir)
		self.rawvcf=rawvcf
		self.chrom=chrom
		self.mpileup=mpileup
		self.software=software
		
	
	def samToolsmPileUpbcfToolsCall(self):
		if self.mpileup:   #for varscan2
			self.mPileup=os.path.join(self.outdir, self.bamBasename.replace('.bam','') + '.samtools.mpileup')
			codeSamtools=[self.software.samtools + ' mpileup',
				'-B ',   #disable BAQ computation. This adjustment seems too stringent, and often reduces base quality scores too aggressively.
				'-q 1 ', #To limit the pileup to reads with mapping quality > 0 (recommended)
				'-f ' + self.software.reference]
			if self.target:
				codeSamtools+=['-l ' + self.target]
			codeSamtools+=[self.bam,
				'>',
				self.mPileup]
			return '   \\\n\t'.join(codeSamtools)
		
		self.rawVcf=os.path.join(self.outdir, self.bamBasename.replace('.bam','') + '.samtools.raw.vcf')
		codeSamtools=[self.software.samtools + ' mpileup',
			'-L 10000 ', #Skip INDEL calling if the average per-sample depth is above INT. [250] for -g or -v
			'-q 30 ', #Minimum mapping quality for an alignment to be used [0]
			'-Q 30 ', #Minimum base quality for a base to be considered [13] 
			'-C 50 ', #Coefficient for downgrading mapping quality for reads containing excessive mismatches. Given a read with a phred-scaled probability q of being generated from the mapped position, the new mapping quality is about sqrt((INT-q)/INT)*INT. A zero value disables this functionality; if enabled, the recommended value for BWA is 50. [0]
			'-t DP,DP4,SP,DV ', #Comma-separated list of FORMAT and INFO tags to output (case-insensitive) with -g or -v
			'-m 2 ', #Minimum number gapped reads for indel candidates INT. [1]  for -g or -v
			'-F 0.002 ', #Minimum fraction of gapped reads [0.002] 
			'-u ', #Generate uncompressed VCF/BCF output, which is preferred for piping
			'-g ',  #Compute genotype likelihoods and output them in the binary call format (BCF). As of v1.0, this is BCF2 which is incompatible with the BCF1 format produced by previous (0.1.x) versions of samtools
			'-f ' + self.software.reference]  #The faidx-indexed reference file in the FASTA format. The file can be optionally compressed by bgzip. [null] 			
		if self.target:
			codeSamtools+=['-l ' + self.target]  #BED or position list file containing a list of regions or sites where pileup or BCF should be generated. If BED, positions are 0-based half-open
		if self.chrom:
			self.rawVcf=os.path.join(self.outdir, self.bamBasename.replace('.bam','') + '.samtools.' + self.chrom + '.raw.vcf')
			codeSamtools+=['-r ' + self.chrom] #Only generate pileup in region. Requires the BAM files to be indexed. If used in conjunction with -l then considers the intersection of the two requests. STR [all sites] 
		codeSamtools+=[self.bam]
		
		codeBcftools=['| '+ self.software.bcftools + ' call', #This command replaces the former bcftools view caller
			'-m ', #alternative modelfor multiallelic and rare-variant calling designed to overcome known limitations in -c calling model (conflicts with -c) 
			'-O v', #Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v)
			'-o ' + self.rawVcf]  #output FILE			
		if not self.rawvcf:
			codeBcftools+=['-v '] #output variant sites only	
		return '   \\\n\t'.join(codeSamtools+codeBcftools)
		

	def vcfFilter(self):
		self.rawfiltVcf=self.rawVcf.replace('.raw.vcf','') + '.rawfilt.vcf'
		self.filtVcf=self.rawfiltVcf.replace('.rawfilt.vcf','') + '.filt.vcf'
		self.snpVcf=self.filtVcf.replace('.vcf','') + '.snp.vcf'
		self.indelVcf=self.filtVcf.replace('.vcf','') + '.indel.vcf'
		bcfFilter=[self.software.bcftools + ' filter',
			'-s FILTER',  #annotate FILTER column with STRING or, with +, a unique filter name generated by the program ("Filter%d")
			'-i "%QUAL>30 && DP>4 && MQ>40"', #include only sites for which EXPRESSION is true
			self.rawVcf,
			'>',
			self.rawfiltVcf]
		passFilter= 'awk -F "\\t" \'{if(/^#/){print}else{if($7=="PASS"){print}}}\' %s > %s ' % (self.rawfiltVcf, self.filtVcf)
		snpFilter= "awk '/^#/ || !/INDEL;/' %s > %s" % (self.filtVcf, self.snpVcf)
		indelFilter= "awk '/^#/ || /INDEL;/' %s > %s" % (self.filtVcf, self.indelVcf)
		return '\n'.join(['   \\\n\t'.join(bcfFilter), passFilter, snpFilter, indelFilter])
		
	
	

	def runJob(self):
		subprocess.check_call(self.samToolsmPileUpbcfToolsCall(), shell=True)
		if not self.mpileup:
			subprocess.check_call(self.vcfFilter(), shell=True)
		
	def __call__(self):
		if self.mpileup:
			return '\n\n'.join([self.samToolsmPileUpbcfToolsCall()]) + '\n'
		else:
			return '\n\n'.join([self.samToolsmPileUpbcfToolsCall(), self.vcfFilter()]) + '\n'


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--bam',help="mapping file",required=True)
	parser.add_argument('--target',help="target region", default=None)
	parser.add_argument('--outdir',help="output path",required=True)
	parser.add_argument('--rawvcf', type=bool, help="whether the middle file rawvcf need to be generated", default=False)
	parser.add_argument('--chrom', help="SNP/INDEL calling for individial chromosome", default=None)
	parser.add_argument('--mpileup', type=bool, help="whether the mpileup file need to be generated", default=False)
	return parser.parse_args()	


	
if __name__ == '__main__':
	args=paramsParse()
	samTools=SAMTOOLS(args.bam, args.outdir, args.target, args.chrom, args.rawvcf, args.mpileup)
	samTools.runJob()