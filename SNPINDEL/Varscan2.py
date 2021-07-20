# coding=utf-8
# independent file to encapsulate the software with nesserray format transfor for use independently
import subprocess
import argparse
import os.path
import sys
import os
import Samtools

usage='''
For somatic SNV/INDEL mutation only

samtools mpileup -B -q 1


#################################
Germline Variant Calling
VarScan parses the pileup input one base at a time, computing the number of bases supporting each observed allele. Only bases meeting the minimum base quality (default: 20) from reads meeting the minimum mapping quality (default: 1) are considered. The coverage (number of qualifying bases) is calculated. If this meets the minimum threshold (default: 20), VarScan examines each allele that was observed, testing to see if it:

    Was supported by at least the minimum number of supporting reads [--min-reads2]
    Meets the minimum allele frequency threshold [--min-var-freq]
    Passes a basic strand-bias filter (if --strand-filter set to 1)
    Has a Fisher's Exact Test p-value below the threshold (if --p-value specified)

For the mpileup2snp/pileup2snp and mpileup2indel/pileup2indel commands, only variants meeting the above criteria are reported. For the mpileup2cns/pileup2cns command, the most-frequent variant meeting the above criteria is reported, or else the reference base is reported. Variants are reported as heterozygotes unless the variant allele frequency exceeds the --min-freq-for-hom which defaults to 0.80. 


filter
This command filters variants in a file by coverage, supporting reads, variant frequency, or average base quality. It is for use with output from pileup2snp or pileup2indel.

	USAGE: java -jar VarScan.jar filter [variants file] OPTIONS
        variants file - A file of SNP or indel calls from VarScan pileup2snp or pileup2indel

	OPTIONS:
	--min-coverage	Minimum read depth at a position to make a call [10]
	--min-reads2	Minimum supporting reads at a position to call variants [2]
	--min-strands2	Minimum # of strands on which variant observed (1 or 2) [1]
	--min-avg-qual	Minimum average base quality for variant-supporting reads [20]
	--min-var-freq	Minimum variant allele frequency threshold [0.20]
	--p-value	Default p-value threshold for calling variants [1e-01]
	--indel-file	File of indels for filtering nearby SNPs, from pileup2indel command
	--output-file	File to contain variants passing filters

##################################


At every position where both normal and tumor have sufficient coverage, a comparison is made. First, normal and tumor are called independently using the germline consensus calling functionality. Then, their genotypes are compared by the following algorithm:

    If tumor matches normal:
        If tumor and normal match the reference
            ==> Call Reference
        Else tumor and normal do not match the reference
            ==> Call Germline

    Else tumor does not match normal:
    Calculate significance of allele frequency difference by Fisher's Exact Test
        If difference is significant (p-value < threshold):
            If normal matches reference
                ==> Call Somatic
            Else If normal is heterozygous
                ==> Call LOH
            Else normal and tumor are variant, but different
                ==> Call IndelFilter or Unknown
        Else difference is not significant:
        Combined tumor and normal read counts for each allele. Recalculate p-value.
            ==> Call Germline 


USAGE: java -jar VarScan.jar somatic [normal_pileup] [tumor_pileup] [output] OPTIONS
OPTIONS:
	--output-snp - Output file for SNP calls [default: output.snp]
	--output-indel - Output file for indel calls [default: output.indel]
	--min-coverage - Minimum coverage in normal and tumor to call variant [8]
	--min-coverage-normal - Minimum coverage in normal to call somatic [8]
	--min-coverage-tumor - Minimum coverage in tumor to call somatic [6]
	--min-var-freq - Minimum variant frequency to call a heterozygote [0.10]
	--min-freq-for-hom	Minimum frequency to call homozygote [0.75]
	--normal-purity - Estimated purity (non-tumor content) of normal sample [1.00]
	--tumor-purity - Estimated purity (tumor content) of tumor sample [1.00]
	--p-value - P-value threshold to call a heterozygote [0.99]
	--somatic-p-value - P-value threshold to call a somatic site [0.05]
	--strand-filter - If set to 1, removes variants with >90% strand bias
	--validation - If set to 1, outputs all compared positions even if non-variant
	
java -jar VarScan.jar processSomatic [output.snp]


The above command will produce 4 output files:

output.snp.Somatic.hc 	(high-confidence Somatic mutations)
output.snp.Somatic.lc 	(low-confidence Somatic mutations)
output.snp.Germline 	(sites called Germline)
output.snp.LOH 		(sites called loss-of-heterozygosity, or LOH)



somaticFilter
This command filters somatic mutation calls to remove clusters of false positives and SNV calls near indels. Note: this is a basic filter. More advanced filtering strategies consider mapping quality, read mismatches, soft-trimming, and other factors when deciding whether or not to filter a variant. See the VarScan 2 publication (Koboldt et al, Genome Research, Feb 2012) for details.

	USAGE: java -jar VarScan.jar somaticFilter [mutations file] OPTIONS
        mutations file - A file of SNVs from VarScan somatic

        OPTIONS:
        --min-coverage  Minimum read depth [10]
        --min-reads2    Minimum supporting reads for a variant [2]
        --min-strands2  Minimum # of strands on which variant observed (1 or 2) [1]
        --min-avg-qual  Minimum average base quality for variant-supporting reads [20]
        --min-var-freq  Minimum variant allele frequency threshold [0.20]
        --p-value       Default p-value threshold for calling variants [1e-01]
        --indel-file    File of indels for filtering nearby SNPs
        --output-file   Optional output file for filtered variants
		


fpfilter

USAGE: java -jar VarScan.jar fpfilter [variant file] [readcount file] OPTIONS
	variant file - A file of SNPs or indels in VarScan-native or VCF format
	readcount file - The output file from bam-readcount for those positions
	***For detailed filtering instructions, please visit http://varscan.sourceforge.net***

	OPTIONS:
	--output-file		Optional output file for filter-pass variants
	--filtered-file		Optional output file for filter-fail variants
	--dream3-settings	If set to 1, optimizes filter parameters based on TCGA-ICGC DREAM-3 SNV Challenge results
	--keep-failures		If set to 1, includes failures in the output file

	FILTERING PARAMETERS:
	--min-var-count		Minimum number of variant-supporting reads [4]
	--min-var-count-lc	Minimum number of variant-supporting reads when depth below somaticPdepth [2]
	--min-var-freq		Minimum variant allele frequency [0.05]
	--max-somatic-p		Maximum somatic p-value [0.05]
	--max-somatic-p-depth	Depth required to test max somatic p-value [10]
	--min-ref-readpos	Minimum average read position of ref-supporting reads [0.1]
	--min-var-readpos	Minimum average read position of var-supporting reads [0.1]
	--min-ref-dist3		Minimum average distance to effective 3' end (ref) [0.1]
	--min-var-dist3		Minimum average distance to effective 3' end (var) [0.1]
	--min-strandedness	Minimum fraction of variant reads from each strand [0.01]
	--min-strand-reads	Minimum allele depth required to perform the strand tests [5]
	--min-ref-basequal	Minimum average base quality for ref allele [15]
	--min-var-basequal	Minimum average base quality for var allele [15]
	--min-ref-avgrl		Minimum average trimmed read length for ref allele [90]
	--min-var-avgrl		Minimum average trimmed read length for var allele [90]
	--max-rl-diff		Maximum average relative read length difference (ref - var) [0.25]
	--max-ref-mmqs		Maximum mismatch quality sum of reference-supporting reads [100]
	--max-var-mmqs		Maximum mismatch quality sum of variant-supporting reads [100]
	--max-mmqs-diff		Maximum average mismatch quality sum (var - ref) [50]
	--min-ref-mapqual	Minimum average mapping quality for ref allele [15]
	--min-var-mapqual	Minimum average mapping quality for var allele [15]
	--max-mapqual-diff	Maximum average mapping quality (ref - var) [50]

'''

class innerSoftware:
	def __init__(self):
		self.java='java'
		self.varscan='VarScan.v2.3.9.jar'
		self.bamreadcount='bam-readcount'
		self.reference='GRCh37.fasta'
		

class VARSCAN:
	def __init__(self, normal, tumor, outdir, purity=1, target=None, software=innerSoftware()):
		self.normal=normal
		self.tumor=tumor
		self.prefix=os.path.basename(self.tumor).replace('.bam','') + '.varscan'
		self.outdir=outdir
		if not os.path.isdir(self.outdir):
			os.mkdir(self.outdir)
		self.purity=str(purity)
		self.software=software
		self.target=target
	
	def formatTrans(self): #from bam to mpileup
		tumorPileup=Samtools.SAMTOOLS(bam=self.tumor, outdir=self.outdir, mpileup=True, target=self.target)	
		normalPileup=Samtools.SAMTOOLS(bam=self.normal, outdir=self.outdir, mpileup=True, target=self.target)
		code=tumorPileup.printJob() + '\n' + normalPileup.printJob()
		self.tumormPileup=tumorPileup.mPileup
		self.normalmPileup=normalPileup.mPileup
		return code
		
#		tumorPileup=Samtools.SAMTOOLS(bam=self.tumor, outdir=self.outdir, mpileup=True, target=self.target)	
#		self.tumormPileup=tumorPileup.mPileup
#		normalPileup=Samtools.SAMTOOLS(bam=self.normal, outdir=self.outdir, mpileup=True, target=self.target)
#		self.normalmPileup=normalPileup.mPileup
#		code=tumorPileup.printJob() + '\n' + normalPileup.printJob()
#		return code
		
		
	def varScan(self):
		code=[self.software.java + ' -jar ' + self.software.varscan + ' somatic',
			self.normalmPileup,
			self.tumormPileup,
			os.path.join(self.outdir, self.prefix),
			'--min-coverage 3',
			'--min-var-freq 0.08',
			'--p-value 0.10',
			'--somatic-p-value 0.05',
			'--tumor-purity ' + self.purity, #Estimated purity (tumor content) of tumor sample [1.00]
			'--normal-purity 1', #Estimated purity (non-tumor content) of normal sample [1.00]
			'--strand-filter 0', #If set to 1, removes variants with >90% strand bias	
			'--output-vcf 1']
			
		self.rawSNP=os.path.join(self.outdir, self.prefix) + '.snp.vcf'
		self.rawIndel=os.path.join(self.outdir, self.prefix) + '.indel.vcf'	
		return '   \\\n\t'.join(code)
		
	def processSomatic(self):
		codeSNP=self.software.java + ' -jar ' + self.software.varscan + ' processSomatic ' + self.rawSNP
		codeIndel=self.software.java + ' -jar ' + self.software.varscan + ' processSomatic ' + self.rawIndel
		self.snpSomatic=os.path.join(self.outdir, self.prefix) + '.snp.Somatic.hc.vcf'
		self.snpGermline=os.path.join(self.outdir, self.prefix) + '.snp.Germline.hc.vcf'
		self.snpLOH=os.path.join(self.outdir, self.prefix) + '.snp.LOH.hc.vcf'
		self.indelSomatic=os.path.join(self.outdir, self.prefix) + '.indel.Somatic.hc.vcf'
		self.indelGermline=os.path.join(self.outdir, self.prefix) + '.indel.Germline.hc.vcf'
		self.indelLOH=os.path.join(self.outdir, self.prefix) + '.indel.LOH.hc.vcf'
				
		return codeSNP + '\n' + codeIndel + '\n'
	
	def somaticFilter(self):
		self.snpSomaticFilter=self.snpSomatic.replace('vcf','') + 'somaticFilter.vcf'
		code=[self.software.java + ' -jar ' + self.software.varscan + ' somaticFilter',
			self.snpSomatic,
			'--min-coverage 10',
			'--min-reads2 2', #Minimum supporting reads for a variant
			'--min-strands2 2',
			'--min-avg-qual 20',
			'--min-var-freq 0.01',
			'--p-value 0.01',
			'--indel-file ' + self.rawIndel,
			'--output-file ' + self.snpSomaticFilter
			]
		
		return '   \\\n\t'.join(code)
	
	
	def _bamReadcount(self, bam, output):   #varscan2 自己也有readcounts的功能VarScan.jar readcounts，可以尝试
		code=[self.software.bamreadcount,
			'-q 1', #minimum mapping quality of reads used for counting
			'-b 15', #minimum base quality at a position to use the read for counting
			bam]
		if self.target:
			code+=[self.target]
		code+=['>', output]
		return ' '.join(code)		
		
	def bamReadcount(self):  #As for bam-readcount, I'd recommend at least mapping quality of 1 and base quality of 15
		self.tumorReadcount=os.path.join(self.outdir, os.path.basename(self.tumor).replace('.bam','') + '.readcount')
		codeTumor=self._bamReadcount(self.tumor, self.tumorReadcount)
		self.normalReadcount=os.path.join(self.outdir, os.path.basename(self.normal).replace('.bam','') + '.readcount')
		codeNormal=self._bamReadcount(self.normal, self.normalReadcount)
		
		return codeTumor + '\n' + codeNormal + '\n'		
	
		
	def _fpFilter(self, varfile, readcount, output):
		code=[self.software.java + ' -jar ' + self.software.varscan + ' fpfilter ',
			varfile,
			readcount,
			'--min-var-count = 3',
			'--min-var-count-lc = 1',
			'--min-strandedness = 0',
			'--min-var-basequal = 30',
			'--min-ref-readpos = 0.20',
			'--min-ref-dist3 = 0.20',
			'--min-var-readpos = 0.15',
			'--min-var-dist3 = 0.15',
			'--max-rl-diff = 0.05',
			'--max-mapqual-diff = 10',
			'--min-ref-mapqual = 20',
			'--min-var-mapqual = 30',
			'--max-var-mmqs = 100',
			'--max-ref-mmqs = 50',
			'--output-file ' + output]
		
		return '   \\\n\t'.join(code)
		
	
	def fpFilter(self):
		self.snpSomaticFilterfpFilter=self.snpSomaticFilter.replace('vcf','') + 'fpFilter.vcf'
		self.snpGermlinefpFilter=self.snpGermline.replace('vcf','') + 'fpFilter.vcf'
		self.snpLOHfpFilter=self.snpLOH.replace('vcf','') + 'fpFilter.vcf'
		
		self.indelSomaticfpFilter=self.indelSomatic.replace('vcf','') + 'fpFilter.vcf'
		self.indelGermlinefpFilter=self.indelGermline.replace('vcf','') + 'fpFilter.vcf'
		self.indelLOHfpFilter=self.indelLOH.replace('vcf','') + 'fpFilter.vcf'
		
		codesnpSomaticFilterfpFilter=self._fpFilter(self.snpSomaticFilter, self.tumorReadcount, self.snpSomaticFilterfpFilter)
		codesnpGermlinefpFilter=self._fpFilter(self.snpGermline, self.normalReadcount, self.snpGermlinefpFilter)
		codesnpLOHfpFilter=self._fpFilter(self.snpLOH, self.normalReadcount, self.snpLOHfpFilter)
		
		codeindelSomaticfpFilter=self._fpFilter(self.indelSomatic, self.tumorReadcount, self.indelSomaticfpFilter)
		codeindelGermlinefpFilter=self._fpFilter(self.indelGermline, self.normalReadcount, self.indelGermlinefpFilter)
		codeindelLOHfpFilter=self._fpFilter(self.indelLOH, self.normalReadcount, self.indelLOHfpFilter)
		
		return '\n\n'.join([codesnpSomaticFilterfpFilter, codesnpGermlinefpFilter, codesnpLOHfpFilter, codeindelSomaticfpFilter, codeindelGermlinefpFilter, codeindelLOHfpFilter]) + '\n'
		
		
	def __call__(self):
		subprocess.check_call(self.printJob(), shell = True)
#		subprocess.check_call(self.formatTrans(), shell=True)
#		subprocess.check_call(self.varScan(), shell=True)
#		subprocess.check_call(self.processSomatic(), shell=True)
#		subprocess.check_call(self.somaticFilter(), shell=True)
#		subprocess.check_call(self.bamReadcount(), shell=True)
#		subprocess.check_call(self.fpFilter(), shell=True)
		
	def printJob(self):
		return '&&\n'.join([self.formatTrans(), self.varScan(), self.processSomatic(), self.somaticFilter(), self.bamReadcount(), self.fpFilter()]) + '\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--normal',help="normal bam file",required=True)
	parser.add_argument('--tumor',help="tumor bam file",required=True)
	parser.add_argument('--outdir',help="path for output",required=True)
	parser.add_argument('--purity',help="purity for tumor sample",default=1, type=int)
	parser.add_argument('--target',help="bed target file",default=None)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	varScan=VARSCAN(args.normal, args.tumor, args.outdir, args.purity, args.target)
	varScan()