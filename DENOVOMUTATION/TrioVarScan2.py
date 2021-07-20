# coding=utf-8
# guoyang
import subprocess
import argparse
import os.path
import sys
import os

usage='''
De novo mutation数量标准
Based on recent WGS studies in families, we think that the de novo mutation rate in humans is approximately 1.1 × 10-8 per haploid genome (1000 Genomes Project Consortium, 2010; Roach et al., 2010). 
By this estimate, an individual's diploid genome harbors, on average, around 64 de novo mutations among 3.2 billion base pairs. 
In the consensus coding sequence (CCDS) (~34 mbp), we expect less than one de novo coding mutation per diploid individual. 


Trio Calling Algorithm
VarScan first calls variants across all three samples the same fashion as it does for mpileup2snp using the default (or user-provided) --min-var-freq and --p-value settings. 
Next, it identifies any variants with apparent Mendelian Inheritance Errors (i.e. present in child but absent from either parent). 
In these instances, it re-calls the samples that should have a variant but were called wild type with adjusted settings (--adj-var-freq and --adj-p-value), 
in an attempt to correct the call. This often corrects the MIE, in which case the corrected genotypes are reported. 
If not, the site will be flagged as mendelError (in the FILTER field) and/or DENOVO (in the INFO field). 


Output
The above command will produce two VCF output files: 
one for SNPs (trio.mpileup.output.snp.vcf) and one for indels (trio.mpileup.output.indel.vcf). 
Relevant INFO fields include:

    FILTER - mendelError if MIE, otherwise PASS
    STATUS - 1=untransmitted, 2=transmitted, 3=denovo, 4=MIE
    DENOVO - if present, indicates a high-confidence de novo mutation call


Downstream Filtering/Interpreation
Even the stringent methodology in VarScan trio will likely call candidate de novo mutations that are not real. 
These should be filtered for false positives (using the child BAM) the same way that somatic mutations are. 
They should also be filtered to aggressively remove known germline variants (i.e. by removing common dbSNPs). 
A recent whole-genome sequencing study in Dutch families developed a random forest classifier for distinguishing true de novo mutations, 
and found that the most important factors were related to sequencing depth and read counts. 
In other words, the ideal de novo mutation call will have high depth (>20x) in all three samples, 
with good support in the child (40-50% of reads for autosomal calls) and no variant-supporting reads in either parent.
'''

class innerSoftware:
	def __init__(self):
		self.java = 'java'
		self.samtools = ''
		self.varscan2 = ''

class innerDatabase:
    def __init__(self):
        self.reference = ''



class TRIOVARSCAN:
    def __init__(self, father, mother, child, outdir, software=innerSoftware(), database = innerDatabase()):
        self.father = os.path.abspath(father)
        self.mother = os.path.abspath(mother)
        self.child = os.path.abspath(child)
        self.outdir = os.path.abspath(outdir)
        self.software = software
        self.database = database
		
    def SamtoolsMpileUp(self):
        self.triompileup = os.path.join(self.outdir, 'trio.mpileup')
        code = [self.software.samtools + ' mpileup',
            '-B',
            '-q 1',
            '-f ' + self.database.reference,
            self.father,
            self.mother,
            self.child,
            '>',
            self.triompileup
        ]		
        return '   \\\n\t'.join(code)

    def TrioVarscan(self):
        self.trioprefix = os.path.join(self.outdir, 'trio.mpileup.output')
        self.trioSNP = self.trioprefix + '.snp.vcf'
        self.trioINDEL = self.trioprefix + '.indel.vcf'
        code = [self.software.java + ' -jar ' + self.software.varscan2 + ' trio',
            self.triompileup,
            self.trioprefix,
            '--min-coverage 10',
            '--min-var-freq 0.20',
            '--p-value 0.05',
            '-adj-var-freq 0.05',
            '-adj-p-value 0.15'
        ]
        return '   \\\n\t'.join(code)
		
    def printJob(self):
        return self.SamtoolsMpileUp() + '\n\n' + self.TrioVarscan()

    def __call__(self):
        subprocess.check_call(self.printJob(), shell=True)
		


def paramsParse():
    parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--father', help = "father bam file", required=True)
    parser.add_argument('--mother', help = "mother bam file", required = True)
    parser.add_argument('--child', help = "child bam file", required = True)
    parser.add_argument('--outdir', help = "output dir", default = '.')
    return parser.parse_args()

if __name__ == '__main__':
    args=paramsParse()
    trioVarscan=TRIOVARSCAN(args.father, args.mother, args.child, args.outdir)
    trioVarscan()