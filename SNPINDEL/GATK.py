# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
the HaplotypeCaller relies on a ploidy assumption (diploid by default) to inform its genotype likelihood and 
variant quality calculations.


sample1_lane1.fq
sample1_lane2.fq 
sample2_lane1.fq 
sample2_lane2.fq

sample1_lane1.dedup.realn.recal.bam 
sample1_lane2.dedup.realn.recal.bam 
sample2_lane1.dedup.realn.recal.bam  
sample2_lane2.dedup.realn.recal.bam 

sample1.merged.bam 
sample2.merged.bam 

sample1.merged.dedup.realn.bam 
sample2.merged.dedup.realn.bam


需要的注释数据库都可以从GATK官网得到
https://software.broadinstitute.org/gatk/download/bundle


Recommended sets of known sites per tool
Summary table
Tool									dbSNP 129		dbSNP >132		Mills indels		1KG indels		HapMap		Omni
RealignerTargetCreator														X					X 		
IndelRealigner																X					X 		
BaseRecalibrator											X				X					X 		
(UnifiedGenotyper/ HaplotypeCaller)							X 				
VariantRecalibrator											X				X									X		X
VariantEval									X 			

##########################################
Variants Calling

java -jar GenomeAnalysisTK.jar \
	-R reference.fasta \
	-T HaplotypeCaller \
	-I sample1.bam \
	--emitRefConfidence GVCF \
	[--dbsnp dbSNP.vcf] \          #加了这个参数，rs号可以出现在ID列
	[-L targets.interval_list] \
	-o output.raw.snps.indels.g.vcf

Joint Genotyping

java -jar GenomeAnalysisTK.jar \
	-T GenotypeGVCFs \
	-R reference.fasta \
	--variant sample1.g.vcf \
	--variant sample2.g.vcf \
	-o output.vcf

#######################################
Single-sample GVCF calling on DNAseq (for `-ERC GVCF` cohort analysis workflow)

   java -jar GenomeAnalysisTK.jar \
     -R reference.fasta \
     -T HaplotypeCaller \
     -I sample1.bam \
     --emitRefConfidence GVCF \
     [--dbsnp dbSNP.vcf] \
     [-L targets.interval_list] \
     -o output.raw.snps.indels.g.vcf
 

 Variant-only calling on DNAseq

   java -jar GenomeAnalysisTK.jar \
     -R reference.fasta \
     -T HaplotypeCaller \
     -I sample1.bam [-I sample2.bam ...] \
     [--dbsnp dbSNP.vcf] \
     [-stand_call_conf 30] \
     [-L targets.interval_list] \
     -o output.raw.snps.indels.vcf
 
 #########################################
Read filters

These Read Filters are automatically applied to the data by the Engine before processing by HaplotypeCaller.

    HCMappingQualityFilter
    MalformedReadFilter
    BadCigarFilter
    UnmappedReadFilter
    NotPrimaryAlignmentFilter
    FailsVendorQualityCheckFilter
    DuplicateReadFilter
    MappingQualityUnavailableFilter

'''

class innerSoftware:
	def __init__(self):
		self.java=''
		self.gatk=''
		
		
class innerDatabase:
	def __init__(self):
		self.reference='human_g1k_v37_decoy.fasta'


class GATK:
	def __init__(self, infile, outdir, regions, software=innerSoftware(), database=innerDatabase()):
		self.infile=infile
		self.prefix=os.path.basename(self.infile).replace('.bam', '')
		self.outdir=outdir
		self.regions=regions
		if self.regions:
			self.suffix＝os.path.basename(self.regions)
		else:
			self.suffix=''
		self.software=software
		self.database=database
		
	def haplotypeCaller(self):
		'''
		The terms "called" and "filtered" are tricky because they can mean different things depending on context. In ordinary language, people often say a site was called if it was emitted as variant. But in the GATK's technical language, saying a site was called means that that site passed the confidence threshold test. For filtered, it's even more confusing, because in ordinary language, when people say that sites were filtered, they usually mean that those sites successfully passed a filtering test. However, in the GATK's technical language, the same phrase (saying that sites were filtered) means that those sites failed the filtering test. In effect, it means that those would be filtered out if the filter was used to actually remove low-confidence calls from the callset, instead of just tagging them. In both cases, both usages are valid depending on the point of view of the person who is reporting the results. So it's always important to check what is the context when interpreting results that include these terms.
		
		Although you now have a nice fresh set of variant calls, the variant discovery stage is not over. The distinctions made by the caller itself between low-confidence calls and the rest is very primitive, and should not be taken as a definitive guide for filtering. The GATK callers are designed to be very lenient in calling variants, so it is extremely important to apply one of the recommended filtering methods (variant recalibration or hard-filtering), in order to move on to downstream analyses with the highest-quality call set possible.
		'''
		self.rawVariants=os.path.join(self.outdir, self.prefix + '.' + self.suffix + '.GATK.raw.vcf')
		code=[self.software.java + ' -Xmx6g -jar ' + self.software.gatk,
			'-T HaplotypeCaller',
			'-R ' + self.database.reference,
			'-I ' + self.infile,
			'-nct 8',
			'--genotyping_mode DISCOVERY',  #This specifies how we want the program to determine the alternate alleles to use for genotyping. In the default DISCOVERY mode, the program will choose the most likely alleles out of those it sees in the data. In GENOTYPE_GIVEN_ALLELES mode, the program will only use the alleles passed in from a VCF file (using the -alleles argument). This is useful if you just want to determine if a sample has a specific genotype of interest and you are not interested in other alleles. 
			'-stand_emit_conf 10',  #This is the minimum confidence threshold (phred-scaled) at which the program should emit sites that appear to be possibly variant.
			'-stand_call_conf 30',  #This is the minimum confidence threshold (phred-scaled) at which the program should emit variant sites as called. If a site's associated genotype has a confidence score lower than the calling threshold, the program will emit the site as filtered and will annotate it as LowQual. This threshold separates high confidence calls from low confidence calls.
			'-o ' + self.rawVariants
			]
		if self.regions:
			code+=['-L ' + self.regions]	
		return '   \\\n\t'.join(code)
		
	def hardFilter(self):
		'''
		(howto) Apply hard filters to a call set
		https://software.broadinstitute.org/gatk/documentation/article.php?id=2806
		'''
		#This creates a VCF file called raw_snps.vcf, containing just the SNPs from the original file of raw variants.
		self.rawSNP=self.rawVariants.replace('.vcf','') + '.snp.vcf'
		code1=[self.software.java + ' -jar ' + self.software.gatk,
			'-T SelectVariants',
			'-R ' + self.database.reference,
			'-V ' + self.rawVariants,
			'-selectType SNP',
			'-o ' + self.rawSNP   
		]
		
		#SNPs matching any of these conditions will be considered bad and filtered out, i.e. marked FILTER in the output VCF file. The program will specify which parameter was chiefly responsible for the exclusion of the SNP using the culprit annotation. SNPs that do not match any of these conditions will be considered good and marked PASS in the output VCF file. 
		self.filtedSNP=self.rawSNP.replace('.raw.snp.vcf','') + '.filted.snp.vcf'
		code2=[self.software.java + ' -jar ' + self.software.gatk,
			'-T VariantFiltration',
			'-R ' + self.database.reference,
			'-V ' + self.rawSNP,
			'--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"',
			#QualByDepth (QD) 2.0
			#This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-reference samples.

			#FisherStrand (FS) 60.0
			#Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.

			#RMSMappingQuality (MQ) 40.0
			#This is the Root Mean Square of the mapping quality of the reads across all samples.

			#MappingQualityRankSumTest (MQRankSum) 12.5
			#This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for mapping qualities (reads with ref bases vs. those with the alternate allele). Note that the mapping quality rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, i.e. this will only be applied to heterozygous calls.

			#ReadPosRankSumTest (ReadPosRankSum) 8.0
			#This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, i.e. this will only be applied to heterozygous calls.
			'--filterName "my_snp_filter"'，  #For SNPs that failed the filter, the variant annotation also includes the name of the filter. 
			'-o ' + self.filtedSNP			
		]
		self.filtedSNPpass=self.filtedSNP.replace('.vcf','') + '.pass.vcf'
		code3=['grep -v "FILTER" ' + self.filtedSNP + ' > ' + self.filtedSNPpass]
		
		self.rawINDEL=self.rawVariants.replace('.vcf','') + '.indel.vcf'
		code4=[self.software.java + ' -jar ' + self.software.gatk,
			'-T SelectVariants',
			'-R ' + self.database.reference,
			'-V ' + self.rawVariants,
			'-selectType INDEL',
			'-o ' + self.rawINDEL  #containing just the Indels from the original file of raw variants			
		]
		
		self.filtedINDEL=self.rawINDEL.replace('.raw.indel.vcf','') + '.filted.indel.vcf'
		code5=[self.software.java + ' -jar ' + self.software.gatk,
			'-T VariantFiltration',
			'-R ' + self.database.reference,
			'-V ' + self.rawINDEL,
			'--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"',
			#QualByDepth (QD) 2.0
			#This is the variant confidence (from the QUAL field) divided by the unfiltered depth of non-reference samples.

			#FisherStrand (FS) 200.0
			#Phred-scaled p-value using Fisher’s Exact Test to detect strand bias (the variation being seen on only the forward or only the reverse strand) in the reads. More bias is indicative of false positive calls.

			#ReadPosRankSumTest (ReadPosRankSum) 20.0
			#This is the u-based z-approximation from the Mann-Whitney Rank Sum Test for the distance from the end of the read for reads with the alternate allele. If the alternate allele is only seen near the ends of reads, this is indicative of error. Note that the read position rank sum test can not be calculated for sites without a mixture of reads showing both the reference and alternate alleles, i.e. this will only be applied to heterozygous calls.
			'--filterName "my_indel_filter"',
			'-o ' + self.filtedINDEL
		]
		
		self.filtedINDELpass=self.filtedINDEL.replace('.vcf','') + '.pass.vcf'
		code6=['grep -v "FILTER" ' + self.filtedINDEL + ' > ' + self.filtedINDELpass]
		
		return '\n\n'.join(map('   \\\n\t'.join,[code1,code2,code3,code4,code5,code6]))
		
	def VQSR(self):
		'''
		GATK variant quality score recalibration (VQSR) intermittently failed to converge on these single sample whole genome inputs
		http://bcb.io/2015/09/17/hg38-validation/
		'''
		pass
		
	def runJob(self):
		subprocess.check_call(self.haplotypeCaller(), shell=True)
		subprocess.check_call(self.hardFilter(),shell=True)
		
	def __call__(self):
		return '\n\n'.join([self.haplotypeCaller(), self.hardFilter()])
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile', help="bam file",required=True)
	parser.add_argument('--outdir', help="output dir",required=True)
	parser.add_argument('--regions', help="regions file to call", dufault=None)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	gatk=GATK(args.infile, args.outdir, args.regions)
	gatk.runJob()