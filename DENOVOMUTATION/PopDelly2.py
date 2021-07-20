# coding=utf-8
# guoyang
import subprocess
import argparse
import os.path
import sys
import os

usage='''
De novo mutation (SV CNV)
理论上可以做群体calling的CNV和SV工具都可以用来做de novo mutation分析，比如：GenomeStrip delly2等等
这里利用delly2的germline calling过程，完成de novo SV的分析
de novo mutation sv原则上只做 DEL DUP
WGS
本脚本实现群体calling

Germline SV calling

    SV calling is done by sample or in small batches to increase SV sensitivity & breakpoint precision

delly call -t DEL -g hg19.fa -o s1.bcf -x hg19.excl sample1.bam

    Merge SV sites into a unified site list

delly merge -t DEL -m 500 -n 1000000 -o del.bcf -b 500 -r 0.5 s1.bcf s2.bcf ... sN.bcf

    Re-genotype merged SV site list across all samples. This can be run in parallel for each sample.

delly call -t DEL -g hg19.fa -v del.bcf -o s1.geno.bcf -x hg19.excl s1.bam

delly call -t DEL -g hg19.fa -v del.bcf -o sN.geno.bcf -x hg19.excl sN.bam

    Merge all re-genotyped samples to get a single VCF/BCF using bcftools merge

bcftools merge -m id -O b -o merged.bcf s1.geno.bcf s2.geno.bcf ... sN.geno.bcf

    Apply the germline SV filter

delly filter -t DEL -f germline -o germline.bcf merged.bcf
'''

class innerSoftware:
    def __init__(self):
        self.delly2='delly_v0.7.3_parallel_linux_x86_64bit'
        self.bcftools='bcftools'

class innerDatabase:
    def __init__(self):
        self.reference = 'human_g1k_v37_decoy.fasta'
        #self.excl= ''
		

class TRIODELLY:
    def __init__(self, input, child, outdir, type, software=innerSoftware(), database=innerDatabase()):
        self.input = os.path.abspath(input)
        self.samplelist = _inputRead(self.input)
        self.child = os.path.abspath(child)
        self.outdir=os.path.abspath(outdir)
        self.type=type
        self.software=software
        self.database=database
	
    def _inputRead(self):
        with open(self.input) as sl:
            samplelist =[]
            for eachLine in sl:
                samplelist.append(eachLine.strip())
        return samplelist

    def _call(self, sampleBam, sampleBcf, vcfgeno = None):
        code = [self.software.delly2 + ' call ',
            '-t ' + self.type,
            '-g ' + self.database.reference,
            #'-q 20',  #min. paired-end mapping quality
            '-o ' + sampleBcf,
            #'-x ' + self.database.excl,   #Exclude telomere and centromere regions
            '--noindels' #no small InDel calling
            #'--indelsize 500', #max. small InDel size
        ]
        if vcfgeno:
            code += ['-v ' + vcfgeno]
        code += [sampleBam]
        return '   \\\n\t'.join(code)

    def firstCall(self):
        code = []
        for eachSample in self.samplelist:
            code += [_call(eachSample, os.path.join(self.outdir, os.path.basename(eachSample) + '.bcf'))]
        return '\n\n'.join(code)


    def firstMerge(self):
        self.mergedBcf = os.path.join(self.outdir, self.type + '.bcf')
        code = [self.software.delly2 + ' merge ',
            '-t' + self.type,
            '-m 500',
            '-n 1000000',
            '-o ' + self.mergeBcf,
            '-b 500',
            '-r 0.5'
        ]
        for eachSample in self.samplelist:
            code += [os.path.join(self.outdir, os.path.basename(eachSample) + '.bcf')]
        return '   \\\n\t'.join(code)

    def secondCall(self):
        code = []
        for eachSample in self.samplelist:
            code += [_call(eachSample, os.path.join(self.outdir, os.path.basename(eachSample) + '.geno.bcf'), self.mergedBcf)]
        return '\n\n'.join(code)

    def secondMerge(self):
        self.reMergedBcf = os.path.join(self.outdir, self.type + '.merged.bcf')
        code = [self.software.bcftools + ' merge',
            '-m id',
            '-O b',
            '-o ' + self.reMergedBcf
        ]
        for eachSample in self.samplelist:
            code += [os.path.join(self.outdir, os.path.basename(eachSample) + '.geno.bcf')]
        return '   \\\n\t'.join(code)

    def filter(self):
        self.filteredBcf = os.path.join(self.outdir, self.type + '.merged.filterd.bcf')
        filterCode = [self.software.delly2+ ' filter '，
            '-t ' + self.type,
            '-f germline',
            '-o ' + self.filteredBcf,
            self.reMergedBcf
        ]
        self.filteredVcf = os.path.join(self.outdir, self.type + '.merged.filterd.vcf')
        transCode = self.software.bcftools + ' view ' + self.filteredBcf + ' > ' + self.filteredVcf	
        return '   \\\n\t'.join(filterCode) + '\n' + transCode
        #下一步可以根据child的样本名，从vcf中筛选parents没有，child有的位点作为de novo mutation
		
    def printJob(self):
        return '\n\n'.join([self.firstCall(), self.firstMerge(), self.secondCall(), self.secondMerge(), self.filter()])

    def __call__(self):
        subprocess.check_call(self.printJob(), shell=True)



def paramsParse():
    parser=argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('--input', help = "bam file list", required = True)
    #parser.add_argument('--child', help = "child bam file name", default = None)  #后续可用来指定child的样本名，自动筛选
    parser.add_argument('--outdir', help="output dir", default='.')
    parser.add_argument('--type', help="SV type: DEL, DUP, INV, TRA, or INS", choices=['DEL', 'DUP', 'INV', 'TRA', 'INS'], required=True)
    return parser.parse_args()

if __name__ == '__main__':
    args=paramsParse()
    triodelly=TRIODELLY(args.input, args.child, args.outdir, args.type)
    triodelly()