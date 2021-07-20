# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R human_g1k_v37_decoy.fasta \
    -L 10:96000000-97000000 \
    -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \
    -I 7156_snippet.bam \
    -o 7156_realignertargetcreator.intervals

java -Xmx8G -Djava.io.tmpdir=/tmp -jar GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R human_g1k_v37_decoy.fasta \
    -targetIntervals 7156_realignertargetcreator.intervals \
    -known INDEL_chr10_1Mb_b37_1000G_phase3_v4_20130502.vcf \ 
    -I 7156_snippet.bam \
    -o 7156_snippet_indelrealigner.bam
	
	
Why do indel realignment?

Local realignment around indels allows us to correct mapping errors made by genome aligners 
and make read alignments more consistent in regions that contain indels.
Genome aligners can only consider each read independently, and the scoring strategies they use to align reads 
relative to the reference limit their ability to align reads well in the presence of indels. 
Depending on the variant event and its relative location within a read, 
the aligner may favor alignments with mismatches or soft-clips instead of opening a gap 
in either the read or the reference sequence. In addition, the aligner's scoring scheme may use arbitrary tie-breaking, 
leading to different, non-parsimonious representations of the event in different reads.
In contrast, local realignment considers all reads spanning a given position. 
This makes it possible to achieve a high-scoring consensus that supports the presence of an indel event. 
It also produces a more parsimonious representation of the data in the region.
This two-step indel realignment process first identifies such regions where alignments may potentially be improved, 
then realigns the reads in these regions using a consensus model that takes all reads in the alignment context together. 


The resulting BAM reduces false positive SNPs and represents indels parsimoniously. First we use RealignerTargetCreator to identify and create a target intervals list (step 1). Then we perform local realignment for the target intervals using IndelRealigner (step 2). 


Changes to the example realigned record:

    MAPQ increases from 60 to 70. The tool increases each realigned record's MAPQ by ten.
    The CIGAR string, now 72M20I55M4S, reflects the realignment containing a 20bp insertion.
    The OC tag retains the original CIGAR string (OC:Z:110M2I22M1D13M4S) and replaces the MD tag that stored the string for mismatching positions.
    The NM tag counts the realigned record's mismatches, and changes from 8 to 24.

Changes to the realigned read's mate record:

    The MC tag updates the mate CIGAR string (to MC:Z:72M20I55M4S).
    The MQ tag updates to the new mapping quality of the mate (to MQ:i:70).
    The UQ tag updates to reflect the new Phred likelihood of the segment, from UQ:i:100 to UQ:i:68.
'''

class innerSoftware:
	def __init__(self):
		self.java='java'
		self.gatk='GenomeAnalysisTK.jar'
		
class innerDatabase:
	def __init__(self):
		self.reference='human_g1k_v37_decoy.fasta'
		self.Mills1KGindell='Mills_and_1000G_gold_standard.indels.b37.vcf'
		


class REALIGNMENT:
	def __init__(self, infile, outdir, software=innerSoftware(), database=innerDatabase()):
		self.infile=infile
		self.prefix=os.path.basename(self.infile).replace('.bam', '')
		self.outdir=os.path.abspath(outdir)
		self.software=software
		self.database=database
		
	def targetCreator(self):
		self.targetIntervals=os.path.join(self.outdir, self.prefix + '.realignertargetcreator.intervals'
		code=[self.software.java + ' -jar ' + self.software.gatk,
			'-T RealignerTargetCreator',
			'-nt 6',
			'-R ' + self.database.reference,
			#'-L 10:96000000-97000000'  #Otherwise, the tool traverses the entire reference genome
			'-known ' + self.database.Mills1KGindell,
			'-I ' ＋ self.infile,  #The tool samples to a target coverage of 1,000 for regions with greater coverage
			'-o ' + self.targetIntervals  #The file is a text-based one-column list with one interval per row in 1-based coordinates. Header and column label are absent. For an interval derived from a known indel, the start position refers to the corresponding known variant. 
		]		
		return '   \\\n\t'.join(code)
		
	def indelRealigner(self):
		self.realnBam=os.path.join(self.outdir, self.prefix + '.realn.bam')
		code=[self.software.java + ' -Xmx8G' + ' -Djava.io.tmpdir=' + os.path.join(self.outdir, 'tmp') + ' -jar' + self.software.gatk,
			'-T IndelRealigner',
			'-R' + self.database.reference,
			'-targetIntervals ' + self.targetIntervals,
			'-known ' + self.database.Mills1KGindell,
			'-I ' + self.infile,
			'-o ' + self.realnBam	
		]
		return '   \\\n\t'.join(code)
		
		
	def __call__(self):
		subprocess.check_call(self.targetCreator(), shell=True)
		subprocess.check_call(self.indelRealigner(), shell=True)
		
	def printJob(self):
		return '\n\n'.join([self.targetCreator(), self.indelRealigner()])
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile',help="bam file",required=True)
	parser.add_argument('--ourdir',help="output dir", default='.')
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	reAlign=REALIGNMENT(args.infile)
	reAlign()