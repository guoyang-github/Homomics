# coding=utf-8

import subprocess
import argparse
import os.path
import sys
import os

usage='''
Control-FREEC is a tool for detection of copy-number changes and allelic imbalances (including LOH) using deep-sequencing data developed by the Bioinformatics Laboratory of Institut Curie (Paris).

Control-FREEC automatically computes, normalizes, segments copy number and beta allele frequency (BAF) profiles, then calls copy number alterations and LOH. 
The control (matched normal) sample is optional for whole genome sequencing data but mandatory for whole exome or targeted sequencing data. 
For whole genome sequencing data analysis, the program can also use mappability data (files created by GEM).

Starting from version 8.0, we provide a possibility to detect subclonal gains and losses and evaluate the likeliest average ploidy of the sample. 
Also, the evaluation procedure for the level of contamination by normal cells has been improved. 




Control-FREEC's output
 _CNVs: file with coordinates of predicted copy number alterations. Information for each region:

    chromosome
    start position
    end position
    predicted copy number
    type of alteration
    genotype (if [BAF] option is used)
    * precentage of uncertainty of the predicted genotyp (max=100).    数字越小越准(The lower the uncertainty, the more reliable is the prediction. The maximal value should be 100. Values lower than 1 are usually trustable)，但－1是不准的(-1 means « information is not available ») 详见ucsc genome browser
    *&** status (somatic or germline)
    *&** precentage of germline CNA/LOH 

 _ratio.txt: file with ratios and predicted copy number alterations for each window. Information for each window:

    chromosome
    start position
    ratio     相对于polidy的比例  or  case/control的比率         区间要么是WGS的windows，要么是wes target区的一个个区间
    median ratio for the whole fragment resulted from segmentation
    predicted copy number
    * median(abs(BAF-0.5)) per window
    * estimated BAF
    * genotype
    * precentage of uncertainty of the predicted genotype, max=100. 

* _BAF.txt: file B-allele frequencies for each possibly heterozygous SNP position. Information for each window:

    chromosome
    start position
    BAF
    fitted frequency of A-allele
    fitted frequency of B-allele
    inferred frequency of A-allele
    inferred frequency of B-allele
    precentage of uncertainty of the predicted genotype, max=100. 

 _sample.cnp and _control.cnp files: files with raw copy number profiles. Information for each window:

    chromosome
    start position
    number of read starts 

 GC_profile.cnp: file with GC-content profile. Information for each window:

    chromosome
    start position
    GC-content
    percentage of ACGT-letter per window (1-poly(N)%)
    *** percentage of uniquely mappable positions per window 

**** _ratio.BedGraph: file with ratios in BedGraph format for visualization in the UCSC genome browser. The file contains tracks for normal copy number, gains and losses, and copy neutral LOH (*). Information for each window:

    chromosome
    start position
    end position
    ratio*ploidy 

*   if [BAF] option is used
**  if [control] option is used
*** if gemMappabilityFile is used
**** if [general] BedGraphOutput=TRUE (see config file options) 





copy number图（数据从ratio文件里来）：
给定一个最大可以显示的Ratio的值maxLevelToPlot
绿色点： ratio$Ratio*ploidy
砖红色： copynumber大于ploidy的点的ratio$Ratio*ploidy，局部的ratio值，整体的copynumber，所以会有低于ploidy的红点儿
中蓝色： copynumber小于ploidy且不为－1的点的ratio$Ratio*ploidy


BAF图（数据从BAF文件里来）：
橘色： A等于0.5的BAF值      (A是文件中的值)
蓝色： A不等于0.5，且不为－1的点的BAF值
黑线： 至少跟后一位A值连续的点的A值，至少跟后一位B值连续的点的B值，部分黑线（0.5处）被紫线盖住了
紫线： 函数拟合的fittedA及fittedB值，连接方式等同于黑线


In general, blue are losses, red are gains, normal copy number is green. 
For the BAF profiles: AB and AA/BB (BAF=0.5) is orange, otherwise it's blue. 
If you uncomment a couple of lines in "makeGraph.R" you could see predicted copy number in black and segment medians in purple


_ratio.BedGraph的内容和_ratio.txt的内容一致
_BAF.txt 都是包含在已知的SNP位点中的
如果config中不使用window及step，而使用coefficientOfVariation, _ratio.txt的区间大小是不一致的
_CNVs中包涵 LOH 及 拷贝中性CN-LOH 看倍性及genotype定


the p-values to predicted CNVs by running assess_significance.R
precentage of uncertainty of the predicted genotype (max=100). (!) It is **not** the uncertainty of a CNA
这是针对CNV和genotype的两个不同指标


WES必须成对样本，WGS可以单个样本
成对样本的情况，并不太需要GC矫正和mapping质量矫正
单个样本的情况，需要chrFiles序列文件计算GC-content profile 或提供GCcontentProfile
单个样本的情况，如果不提供gemMappabilityFile，就用N fraction
单个样本也可以得到LOH（这里应该还要考虑 1.本样本自己的GATK的结果 还是 2.dbSNP里的位点结果的不同）


(For CNV) The algorithm includes several steps. 
First, it calculates the raw CNP by counting reads in non-overlapping windows. 
If not provided by the user, window size can be automatically selected 
using depth of coverage infor- mation to optimize accuracy of CNA prediction. 

The second step is profile normalization. 
If a control is not provided by the user, we compute the GC- content profile. 
The normalization procedure of RC by GC-content (or by control RC) is described below. 
The third step is segmentation of the nor- malized CNP. 
To do this we implemented a LASSO-based algorithm sug- gested by (Harchaoui and Lévy-Leduc, 2008). 
Segmentation provided by this algorithm is robust against outliers, 
which makes it suitable for seg- mentation of deep sequencing CNPs. 
The last step involves analysis of segmented profiles. 
This includes identification of regions of genomic gains 
and losses and prediction of copy number changes in these regions. 
To normalize a raw CNP we fit the observed RC by the GC-content (or the control RC if it is available). 
We base our fitted model on several assumptions: 
(i) the sample main ploidy P is provided, 
(ii) the observed RC in P- copy regions (i.e., regions with copy number equal to P) can be modeled as a polynomial of GC-content (or of control RC), 
(iii) the observed RC in a region with altered copy number is linearly proportional to the RC in P- copy regions, and 
(iv) the interval of measured GC-contents (resp. control RCs when a control dataset is available) in the main ploidy regions must include the interval of all measured GC-contents (resp. control RC). 
The polynomial’s degree is a user defined parameter with a default value of three. 
We provide an initial estimate of the polynomial's parameters 
and then optimize these parameters by iteratively selecting data points 
related to P-copy regions and making a least-square fit on these points only (See Suppl. Methods for more details). 
The resulting polynomial is then used to normalize the CNP (Fig. 1). 
The user has an option to include mappability information into the normalization procedure.


Workflow. The workflow of Control-FREEC consists of 3 steps: 
(1) calculation and segmentation of copy number profiles; 
(2) calculation and segmentation of smoothed B-allele frequency (BAFs) profiles; 
(3) prediction of final genotype status, i.e., copy number and allelic content for each segment (for example, A, AB, AAB, etc.).

'''

class innerSoftware:
	def __init__(self):
		self.controlfreec='/controlFreeC/freec'
		self.makeGraph='/controlFreeC/makeGraph.R'
		self.assess_significance='/controlFreeC/assess_significance.R'
		self.samtools='/samtools/samtools-0.1.18/samtools'

class innerDatabase:
	def __init__(self):
		self.chrFiles='human_b37/byChr/'
		self.chrLenFile='/human_b37/freec/chr.24.length' #可以用reference的.fai文件
		#self.gemMappabilityFile=''
		self.reference='/human_b37/human_g1k_v37_decoy.fasta'
		self.SNPfile='/freec/b37_snp137.SingleDiNucl.1based.txt'
		self.SNPpositions='/hg19_snp142.SingleDiNucl.1based.bed'  #这两个参数可以合并用一个vcf文件了，选用这个样本自己的GATK germline SNP vcf文件
		

class CONTROLFREEC:
	def __init__(self, case, control, outdir, loh, ploidy, sex, nt, type, target, snp = None, software=innerSoftware(), database=innerDatabase()):
		self.case=os.path.abspath(case)
		if control:
			self.control=os.path.abspath(control)
		else:
			self.control=control
		self.outdir=os.path.abspath(outdir)
		self.type=type
		self.target=target
		if self.type == 'WES':
			assert self.control
			assert self.target
		self.loh=loh
		self.ploidy=str(ploidy)
		self.sex=sex
		self.nt=nt
		#self.purity=purity
		#self.contamination=contamination
		#self.format=format
		self.snp=snp
		self.software=software
		self.database=database
		if self.snp:
			self.database.SNPfile=self.snp
			self.database.SNPpositions=self.snp

		
	def controlFreeCConfig(self):
		config=[]
		
		config.append('\n[general]')
		config.append('BedGraphOutput=TRUE')  # if you want an additional output in BedGraph format for the UCSC genome browser
		#bedtools
		#breakPointThreshold      #Default: 0.8   use something like 0.6 to get more segments (and thus more predicted CNVs) 
		config.append('breakPointType = 4')  #desired behavior in the ambiguous regions (poly-N or low mappability regions between two different copy number values) Default: 2
			#0: the "unknown" region is attached to the "known" region on the right
			#1: make a separate fragment of this “unknown” region and then attache it to the left or to the right region choosing the longer one
			#2: make a separate fragment of this “unknown” region and then attache it to the left or to the right region but the “ploidy” copy number has a priority
			#3: make a separate fragment of this “unknown” region and then attache it to the left or to the right region choosing the longer one but this “known” region should make at least half-size of the “unknown” region
			#4: make a separate fragment of this “unknown” region and do not assign any copy number to this region at all 		
		config.append('chrLenFile = ' + self.database.chrLenFile)  #file with chromosome lengths chromosomes that are not in this list won't be considered by Control-FREEC!
		#config.append('coefficientOfVariation = 0.05')  #coefficient of variation to evaluate necessary window size
			#Parameter « coefficientOfVariation » is used only when you don’t provide window size. To use « coefficientOfVariation », comment or delete « window=... » from your config and you will get window size evaluated by Control-FREEC. The default value of « coefficientOfVariation » is 0.05. However, I often use something like « coefficientOfVariation=0.062 » to get smaller windows. 
		config.append('contaminationAdjustment = TRUE')  #set TRUE to correct for contamination by normal cells. If "contamination" is not provided, it will automatically evaluate the level of contamination  Default: FALSE
		#contamination
		#degree
		#intercept
		#minCNAlength
		#minExpectedGC
		#maxExpectedGC
		#minimalSubclonePresence
		if self.type == 'WGS':
			config.append('window=50000')
			#config.append('forceGCcontentNormalization = 0')  #set 1 or 2 to correct the Read Count (RC) for GC-content bias and low mappability even when you have a control sample      Default (WGS): 0  Default (WES): 1 (≥ v9.5) and 0 (< v9.5)
				#0: simply model "sample RC ~ Control RC"
				#1: normalize the sample and the control RC using GC-content and then calculate the ratio "Sample RC/contol RC"
				#2: model "sample RC ~ Control RC" bias, and then normalize for GC-content 
				##with targeted sequencing, I would not recommend to use forceGCcontentNormalization=1 or 2 since capture bias can be much stronger than GC-content bias
			if not self.control: #单样本的时候再对GC及mappability做限制
				#degree   #Default: 3&4 (GC-content based normalization, WGS) or 1 (control-read-count-based normalization, WES)    综合forceGCcontentNormalization这个参数的注释一起看，才能理解freec对GC校正和controlRC校正的操作方式，对GC用多项式的方式（类同于lowess）对controlRC就是算ratio（1阶）
				config.append('chrFiles = ' + self.database.chrFiles)  #path to the directory with chromosomes fasta files (necessary to calculate a GC-content profile if a control dataset and GC-content profile are not available)
				#GCcontentProfile   #GC-content profile for a given window-size (higher priority than chrFiles) Optional! Necessary only if both a control dataset or chromosome sequences (.fasta for the genome of interest) are not available
				#config.append('gemMappabilityFile = ' + self.database.gemMappabilityFile)  #.gem file with information about mappable positions (GEM output)  At this point, gemMappabilityFile can be used only in the mode without control sample. 
				config.append('uniqueMatch = FALSE')  #Use a mappability profile to correct read counts (in this case a mappability file must be provided with "gemMappabilityFile" )   Default: FALSE 
				#minMappabilityPerWindow   #only windows with fraction of mappable positions higher than or equal to this threshold will be considered (if "gemMappabilityFile" is not provided, one uses the percentage of non-N letters per window)  Default: 0.85 
					#If the parameter gemMappabilityFile is not specified, then the fraction of non-N nucleotides per window is used as Mappability. 		
				
		if self.type == 'WES':  #必须是成对样本
			#degree
			config.append('window=0') #Versions < v8.0: We optimized FREEC for window size "window=500" and "step=250"   Versions > v8.0: Use "window=0" to calculate read count per exon
			config.append('chrFiles = ' + self.database.chrFiles)  #根据forceGCcontentNormalization对于WES的默认参数，需要GC校正，所以提供这个参数
			config.append('noisyData = TRUE')  #set TRUE for target resequencing data (e.g., exome-seq) to avoid false positive predictions due to nonuniform capture Default: FALSE 
			config.append('printNA = FALSE')  #set FALSE to avoid printing "-1" to the _ratio.txt files  Useful for exome-seq or targeted sequencing data  Default: TRUE
			config.append('readCountThreshold = 50')  #threshold on the minimal number of reads per window in the control sample  Useful for exome-seq or targeted sequencing data  Default: 10  recommended value >=50 for for exome data 
		config.append('maxThreads = 8')  #number of threads (multi-threading mode)
		config.append('outputDir = ' + self.outdir)
		config.append('ploidy = ' + self.ploidy) #genome ploidy; In case of doubt, you can set different values and Control-FREEC will select the one that explains most observed CNAs
		#sambamba
		#SambambaThreads
		config.append('samtools = ' + self.software.samtools)
		#step
		#telocentromeric
			#The parameters chrLenFile and ploidy are required. Either coefficientOfVariation or window must be specified. Either chrFiles or GCcontentProfile must be specified too if no control dataset is available. 
		if self.sex:
			config.append('sex = ' + self.sex)  #sample sex   "sex=XX" will exclude chr Y from the analysis   "sex=XY" will not annotate one copy of chr X and Y as a loss.
		
		
		config.append('\n[sample]')
		config.append('mateFile = ' + self.case)
		config.append('inputFormat = BAM')
		config.append('mateOrientation = FR')  #0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs) 
		
		if self.control:
			config.append('\n[control]')
			config.append('mateFile = ' + self.control)
			config.append('inputFormat = BAM')
			config.append('mateOrientation = FR')  #0 (for single ends), RF (Illumina mate-pairs), FR (Illumina paired-ends), FF (SOLiD mate-pairs) 

		if self.loh:   #We characterize the allelic content via the B allele frequency(BAF). We limit the list of genomic positions which we consider to evaluate allelic content to known SNPs only.By the B allele, we mean the alternative variant in dbSNP.
			config.append('\n[BAF]')
			config.append('SNPfile = ' +  self.database.SNPfile)  #file with known SNPs
			config.append('minimalCoveragePerPosition = 5')  #minimal read coverage for a position to be considered in BAF analysis
			config.append('minimalQualityPerPosition = 5')  #minimal sequencing quality for a position to be considered in BAF analysis
			config.append('makePileup = ' + self.database.SNPpositions)  #path to a BED/VCF file with SNP positions to create a mini pileup file from the initial BAM file; If provided, a BAM file can be used for the calculation of BAF
				#Optional; if not used, a pileup file should be provided in the [sample] and [control] groups of parameters) This file can be created from hg19_snp142.SingleDiNucl.1based.txt.gz using a command line  gunzip -c hg19_snp142.SingleDiNucl.1based.txt.gz | awk {'printf ("%s\t%s\t%s\n", $1,$2-1,$2)'} >hg19_snp142.SingleDiNucl.1based.bed
			config.append('fastaFile = ' + self.database.reference)  #one fasta file for the whole genome	Optional; need to provide it only when option "makePileup" is used to create a minipileup file from .BAM
			config.append('shiftInQuality = 33')  #basis for Phred quality
			
		if self.target:
			config.append('\n[target]')
			config.append('captureRegions = ' + self.target)  #file with capture regions in .BED format  sorted .BED file should contain the following colomns: chr   0-based start   1-based end

		self.configFile=os.path.join(self.outdir, 'controlFreeC.config')
		with open(self.configFile, 'w') as controlFreeCConfigFile:
			controlFreeCConfigFile.writelines(map(lambda x: x + '\n', config))


	def controlFreeC(self):
		self.CNVs=os.path.join(self.outdir, os.path.basename(self.case) + '_CNVs')
		self.ratio=os.path.join(self.outdir, os.path.basename(self.case) + '_ratio.txt')
		self.BAF=os.path.join(self.outdir, os.path.basename(self.case) + '_BAF.txt')
		code=self.software.controlfreec + ' -conf ' + self.configFile
		return code
		
	def assessSignificance(self):
		self.CNVsPV=self.CNVs + '.p.value.txt'
		code='cat ' + self.software.assess_significance + ' | R --slave --args ' + self.CNVs + ' ' + self.ratio 
		return code
		
	def makeGraph(self):
		self.ratioGraph=self.ratio + '.png'
		code='cat ' + self.software.makeGraph + ' | R --slave --args ' + self.ploidy + ' ' + self.ratio
		if self.loh:
			self.BAFGraph=self.BAF + '.png'
			code += ' ' + self.BAF
		return code
		
	def runJob(self):
		self.controlFreeCConfig()
		subprocess.check_call(self.controlFreeC(), shell=True)
		subprocess.check_call(self.assessSignificance(), shell=True)
		subprocess.check_call(self.makeGraph(), shell=True)
		
	def __call__(self):
		code=[]
		self.controlFreeCConfig()
		code.append(self.controlFreeC())
		code.append(self.assessSignificance())
		code.append(self.makeGraph())
		return '\n\n'.join(code)
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--loh',help="If the LOH included", action='store_true')
	parser.add_argument('--outdir',help="The directory to store software output, default: cwd",default = '.')
	parser.add_argument('--ploidy', help="ploidy", default='2,3,4,5,6')
	parser.add_argument('--sex',help="gender of the sample", choices=['XX','XY'], default=None)
	parser.add_argument('--nt', help="thread number default: 8", type=int, default=6)
	parser.add_argument('--case',help="sample bam file", required=True)
	parser.add_argument('--control',help="control bam file", default=None)
	#parser.add_argument('--purity',help="purity from eg. Absolute", type=float, default=None)
	parser.add_argument('--type',help="sequencing type", required=True, choices=['WGS','WES'])
	parser.add_argument('--target',help="target region for capture sequencing", default=None)
	#parser.add_argument('--contamination',help="going to evaluate contamination", type=bool, choices=[False, True], default=True)
	#parser.add_argument('--format',help="input alignment file format", choices=['BAM','pileup'], default='BAM')
	parser.add_argument('--snp',help="germline snp vcf from GATK", default=None)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	controlFreeC=CONTROLFREEC(args.case, args.control, args.outdir, args.loh, args.ploidy, args.sex, args.nt, args.type, args.target, args.snp)
	controlFreeC.runJob()
	
	
