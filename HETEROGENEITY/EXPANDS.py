# coding=utf-8
# guoyang
# 20170706 完成py封装和R函数   未测试
import subprocess
import argparse
import os.path
import sys
import os

usage='''
EXPANDS predicts tumor purity and subclonal composition from sequencing data.
characterizes genetically diverse subpopulations (SPs) in a tumor using copy number and allele frequencies derived from exome- or whole genome sequencing input data. 

The model predicts subpopulations based on Two assumptions:
1. Two independent driver-events of the same type will not happen at the exact same genomic position in two di↵erent cells. Therefore, no more than two distinct cell populations co-exist with respect to a specific locus.
2. Multiple passenger mutations accumulate in a cell before a driver mutation causes a clonal expansion. Thus, each clonal expansion is marked by multiple mutations.
These two assumptions are translated into the ExPANdS model in five main steps: 
1. cell frequency estimation, 
2. clustering, 
3. filtering, 
4. assignment of mutations to clusters，
5. phylogenetic tree estimation. 
You can demonstrate each of these steps separately. 
The main function runExPANdS performs all five steps. 
The robustness of the subpopulation predictions by ExPANdS increases with the number of mutations provided. It is recommended that at least 200 mutations are used as an input to obtain stable results.
可以分步完成，见《Description of expands》文档，可以一起完成见函数 runExPANdS

#runExPANdS(SNV, CBS, maxScore=2, max_PM=6, min_CellFreq=0.1, precision=NA, plotF=2,snvF=NULL,maxN=8000,region=NA)


对输入数据的要求： somatic位点数：The robustness of the subpopulation predictions by ExPANdS increases with the number of mutations provided. 
It is recommended that at least 200 mutations are used as an input to obtain stable results.
The remaining set of somatic point mutations can be extended to contain loss of heterozygosity (LOH), that is loci with heterozygous germline polymorphisms where the mutated allele is overrepresented in the cancer cell. 
For tumor types with a low number of somatic point mutations, this approach can provide a sufficient number of somatic events for the subsequent procedure.
somatic不够 LOH来凑

对ploidy的理解，这个软件并不是指整套染色体的倍型，还是和CNV类似的概念
主要就是做亚克隆结构的，还能画进化树，纯度就单纯指最大亚克隆的占比
在linux下输出图像还是有点儿问题
'''

class innerSoftware:
	def __init__(self):
        self.python = ''
		self.R = ''
		self.expands = os.path.abspath('./EXPANDS.R')
		

class EXPANDS:
	def __init__(self, somatic, loh, cnv, target, outdir, software=innerSoftware()):
		self.somatic = os.path.abspath(somatic)
		self.loh = loh
		self.cnv = os.path.abspath(cnv)
		self.target = os.path.abspath(target)
		self.outdir = outdir
		self.software = software
		
	def dataTransfor(self):
		self.outSNV = os.path.join(self.outdir, 'varscan2.expands.snv')
		self.outCNV = os.path.join(self.outdir, 'controlfreec.expands.cbs')
		codeSNV = [self.software.python + ' ' + self.software.transformer + ' varscan4EXPANDS',
			'--somatic ' + self.somatic,
			'--loh ' + self.loh,
			'--outdir ' + self.outdir
		]

		codeCNV += [self.software.python + ' ' + self.software.transformer + ' controlfreec4EXPANDS',
			'--cnv ' + self.cnv,
			'--outdir ' + self.outdir
		]
		return '   \\\n\t'.join(codeSNV) + '\n' + '   \\\n\t'.join(codeCNV)
		
	def exPands(self):
		code = [self.software.R + self.software.expands,
			'--snv ' + self.outSNV,
			'--cnv' + self.outCNV,
			'--region' +  self.target,
			'--outdir' + self.outdir
		]
		return '   \\\n\t'.join(code)



	def printJob(self):
		return '&&\n'.join(self.dataTransfor(), self.exPands())

	def __call__(self):
		subprocess.check_call(self.printJob(), shell=True)
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--somatic', help="snpSomaticFilterfpFilter", required=True)
	parser.add_argument('--loh', help="snpLOHfpFilter", default = None)
	parser.add_argument('--cnv', help="controlfreec _CNV", required=True)
	parser.add_argument('--target', help = "coding regions", required = True)
	parser.add_argument('--outdir', help="output dir", default = '.')
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	expands=EXPANDS(args.somatic, args.loh, args.cnv, args.target, args.outdir)
	expands()