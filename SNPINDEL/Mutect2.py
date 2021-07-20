# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine 
of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

MT和HC的区别：
While the HaplotypeCaller relies on a ploidy assumption (diploid by default) to inform its genotype likelihood and 
variant quality calculations, MuTect2 allows for a varying allelic fraction for each variant, 
as is often seen in tumors with purity less than 100%, multiple subclones, and/or copy number variation (either local or aneuploidy). 
MuTect2 also differs from the HaplotypeCaller in that it does apply some hard filters to variants before producing output. 
Finally, some of the parameters used for ActiveRegion determination and graph assembly are set to different default values in MuTect2 compared to HaplotypeCaller.


Tumor/Normal variant calling

   java -jar GenomeAnalysisTK.jar \
     -T MuTect2 \
     -R reference.fasta \
     -I:tumor tumor.bam \
     -I:normal normal.bam \
     [--dbsnp dbSNP.vcf] \
     [--cosmic COSMIC.vcf] \
     [-L targets.interval_list] \
     -o output.vcf
 

Normal-only calling for panel of normals creation

   java -jar GenomeAnalysisTK.jar \
     -T MuTect2 \
     -R reference.fasta \
     -I:tumor normal1.bam \
     [--dbsnp dbSNP.vcf] \
     [--cosmic COSMIC.vcf] \
     --artifact_detection_mode \
     [-L targets.interval_list] \
     -o output.normal1.vcf


For full PON creation, call each of your normals separately in artifact detection mode as shown above. 
Then use CombineVariants to output only sites where a variant was seen in at least two samples:

   java -jar GenomeAnalysisTK.jar \
     -T CombineVariants \
     -R reference.fasta \
     -V output.normal1.vcf -V output.normal2.vcf [-V output.normal2.vcf ...] \
     -minN 2 \
     --setKey "null" \
     --filteredAreUncalled \
     --filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
     [-L targets.interval_list] \
     -o MuTect2_PON.vcf
'''

class innerSoftware:
	def __init__(self):
		self.java = ''
		self.gatk = ''

class innerDatabase:
	def __init__(self):
		self.dbsnp = ''
		self.cosmic = ''
		self.reference = ''

class MUTECT2:
	def __init__(self, tumor, normal, target, output, software=innerSoftware(), database = innerDatabase()):
		self.tumor = tumor
		self.normal = normal
		self.target = target
		self.output = output
		self.software = software
		self.database = database
		
	def MuTect2(self):
		code=[self.software.java + ' -jar ' + self.software.gatk,
			'-T MuTect2',
			'-R ' + self.reference,
			'-I:tumor ' + self.tumor,
			'-I:normal ' + self.normal,
			'--dbsnp ' + self.database.dbsnp,
			'--cosmic ' + self.database.cosmic,
			'-o ' + self.output
		]
		if self.target:
			code += ['-L ' + self.target]
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		subprocess.check_call(self.MuTect2, shell=True)
		
	def printJob(self):
		return self.MuTect2() + '\n\n'
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--tumor', help = "tumor.bam", required = True)
	parser.add_argument('--normal', help = "normal.bam", required = True)
	parser.add_argument('--target', help = "targets.interval_list", default = None)
	parser.add_argument('--output', help = "output.vcf", default = 'output.vcf')
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	muTect2=MUTECT2(args.tumor, args.normal, args.target, args.output)
	muTect2()
