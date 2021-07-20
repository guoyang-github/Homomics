# coding=utf-8
# guoyang
import subprocess
import argparse
import os.path
import sys
import os

usage='''
应用：WES
纯合子区域即基因组上的纯和等位基因的区域。这个区间的产生是由于父母遗传至子代
纯合子定位分析一般用于近亲结婚引起的遗传病，家系间隐性遗传病或认为是 denovo 引起的疾病使用此方法也是有可能得到结果。
对于各类软件的评估，有文献肯定了H3M2在短ROH(<1.5M)上比Plink等有优势，长ROH(>1.5M) 上相差不大。
关于ROH的长度，长度是没有限定的(ranging in length from tens of kilobases to megabases)，其分布也不均匀。


定义：
Runs of Homozygosity (ROH) are chromosomal stretches that in a diploid genome appear in the homozygous state, 
that is, display identical alleles at multiple contiguous loci.


原理及优势：
Currently available sliding window methods, such as PLINK (Purcell et al., 2007) or GERMLINE (Gusev et al., 2009), 
were developed for the analysis of uniformly spaced SNP-array maps and do not fit easily to 
the analysis of the sparse and non uniform distribution of WES SNP maps. 
The identification of long ROH, typically those detected in a context of consanguinity, 
may be faced by applying these methods to WES data (Pippucci et al., 2011). 
However, medium or short ROH cannot be as easily captured from WES data by using traditional sliding windows approaches. 
To meet the need of an approach specifically tailored to WES data, that could overcome the inherent limitations of currently available tools, 
we developed a novel computational approach (homozygosity heterogeneous hidden markov model, H3M2 ) for the identification of ROH. 
The algorithm is based on an heterogeneous hidden markov model that, incorporating the distance between consecutive polymorphic positions into the transition probabilities matrix, 
is able to detect with high sensitivity and specificity ROH of every genomic size. 
The key feature of the algorithm is its heterogeneity, which makes it well-suited for WES data.


科普：
long ROH (IBD IBS)
short ROH
medium-sized ROH

Recent parental relatedness favors the formation of long ROH (several megabases) that occur due to 
IBD (Identity by Descent, when the two alleles at a locus match because they originate from the same common ancestor), 
as opposed to IBS (Identity by State, when the two alleles at a locus match simply by coincidence). 
Homozygosity originating from the occurrence of individual IBD regions due to parental relatedeness (autozygosity) 
is known to possibly contain recessive highly penetrant, deleterious disease-causing mutations surrounded 
by an unusually long homozygous haplotype. 
This is the principle that inspired homozygosity mapping in the study of rare recessive disorders affecting inbred individuals (Lander and Botstein, 1987). 
In outbred individuals, short (up to few hundreds of kilobases) or medium-sized ROH (from hundreds of kilobases to a few megabases) can surround disease-causing mutations as well (Hildebrandt et al., 2009), 
playing a role in predisposition to disease through the effect of mildly deleterious recessive variants (Wang et al., 2009).



使用：
The H3M2 package is made of two modules: 
    The bam parsing module (H3M2BamParsing.sh) that creates the BAF values and 
    the analysis module (H3M2Analyze.sh) that identifies the ROH. 

The first analysis step consists in caluclating the BAF values by using the following command:

/.../H3M2Tool/H3M2BamParsing.sh path2prog path2bamfold filebam path2results ProjectName SampleName hgref FileBed

The second step consists in analyzing the BAF values for detecting ROH by using the following command:

/.../H3M2Tool/H3M2Analyze.sh path2prog path2Results ProjectName SampleName FileBed DNorm P1 P2 Factor

where:

path2prog -> the path to the H3M32 folder
path2bamfold -> the path to the folder with the bam file
filebam -> the name of the bam file
path2results -> the path to the Results folder
ProjectName -> the name of the project (a folder with this label will be created under path2Results)
SampleName -> the name of the sample (a folder with this name will be created under path2Results/projectName)
hgref -> the full path to the human reference genome in .fa (.fasta) format
FileBed -> the full path to the .bed file that contain the snp position to be used for the analysis.
DNorm -> parameter of the H3M2 (it can be 1000, 10000 and 100000. We suggest to use 100000)
P1 -> parameter of the H3M2 (set to 0.1)
P2 -> parameter of the H3M2 (set to 0.1)
Factor -> parameter of the H3M2 (set to 5)


配置文件：
Please note that the H3M2 package contains two snp position .bed files (SNP1000GP.HGb37_Exome.bed and SNP1000GP.HGb37_Exome.mod.bed) 
each one containing the positions of more than 4 milions markers (the markers cover the coding regions of the genome plus flanking regions). 
Use SNP1000GP.HGb37_Exome.mod.bed if the reads have been mapped to the UCSC reference genome (hg19, where chromosomes are encoded as chr1, chr2,....). 
Use SNP1000GP.HGb37_Exome.bed if the reads have been mapped to the 1KGP reference genome (where chromosomes are encoded as 1, 2,....).


结果文件：
The results of the H3M2 analysis will be stored in "path2results/ProjectName/SampleName/Results". 
The folder will contain the ROH identified by H3M2 in .bed format 
(Chromosome, Start, End, Probability of Autozygosity, # of SNPs within the ROH) and a .pdf file 
with the ROH plot for each chromosome (the ROH are highlighted in red ).

'''

class innerSoftware:
	def __init__(self):
		self.path2prog='java'

class innerDatabase:
    def __init__(self):
        self.reference = ''
        self.FileBed = ''
		

class H3M2:
    def __init__(self, infile, outdir, projectname, samplename, software=innerSoftware(), database = innerDatabase()):
        self.inpath = os.path.absdir(infile)
        self.infile = os.path.absname(infile)
        self.outdir = os.path.abspath(outdir)
        self.projectname = projectname
        self.samplename = samplename
        self.software = software
        self.database = database
		
    def H3M2BamParsing(self):
        code = [os.path.join(self.software.path2prog, 'H3M2BamParsing.sh'),
            self.software.path2prog,
            self.inpath,
            self.infile,
            self.outdir,
            self.projectname,
            self.samplename,
            self.database.reference,
            self.database.FileBed
		]		
        return '\t'.join(code) + '\n'

    def H3M2Analyze(self):
        code = [os.path.join(self.software.path2prog, 'H3M2Analyze.sh'),
            self.software.path2prog,
            self.outdir,
            self.projectname,
            self.samplename,
            self.database.FileBed,
            100000,
            0.1,
            0.1,
            5
        ]
        return '\t'.join(code) + '\n'
		
    def printJob(self):
        return self.H3M2BamParsing() + self.H3M2Analyze() + '\n'

    def __call__(self):
        subprocess.check_call(self.printJob(), shell=True)
		


def paramsParse():
    parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--infile', help="bam file path",required=True)
    parser.add_argument('--outdir', help = "output dir", default = '.')
    parser.add_argument('--projectname', help = "project name", default = 'Project')
    parser.add_argument('--samplename', help = "sample name", required = True)
    return parser.parse_args()

if __name__ == '__main__':
    args=paramsParse()
    h3m2=H3M2(args.infile, args.outdir, args.projectname, args.samplename)
    h3m2()