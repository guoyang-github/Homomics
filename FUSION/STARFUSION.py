# coding=utf-8
# guoyang
import subprocess
import argparse
import os.path
import sys
import os

usage='''
https://github.com/STAR-Fusion/STAR-Fusion/wiki
可以从fastq开始，也可以从BAM开始，但是对star的比对参数有特殊设置，所以还要重新比对，不如直接从fastq开始


已经建设好的genome lib 包括 人 和 小鼠
Building a Custom FusionFilter Dataset
https://github.com/FusionFilter/FusionFilter/wiki/Building-a-Custom-FusionFilter-Dataset
需要：
1. genome in FASTA format
2. gene annotation set in GTF format （要包括transcript_id, gene_id, and gene_symbol）



 STAR --genomeDir ${star_index_dir} \                                                                                     
      --readFilesIn ${left_fq_filename} ${right_fq_filename} \                                                                      
      --twopassMode Basic \                                                                                                      
      --outReadsUnmapped None \                                                                                                  
      --chimSegmentMin 12 \                                                                                                    
      --chimJunctionOverhangMin 12 \                                                                                           
      --alignSJDBoverhangMin 10 \                                                                                              
      --alignMatesGapMax 100000 \                                                                                             
      --alignIntronMax 100000 \                                                                                                
      --chimSegmentReadGapMax parameter 3 \                                                                                    
      --alignSJstitchMismatchNmax 5 -1 5 5 \
      --runThreadN ${THREAD_COUNT} \                                                                                                           
      --limitBAMsortRAM 31532137230 \                                                                                           
      --outSAMtype BAM SortedByCoordinate 



## De novo reconstruct fusion transcripts using Trinity
The reconstructed transcripts described above by the 'fusion coding effect' module are based on the reference annotations 
and the reference genome sequence. If you're interested in doing de novo reconstruction of fusion transcripts based 
on the actual fusion-supporting RNA-Seq reads, capturing any additional variant information or novel sequence features 
that may be evidence among the reads, you can include de novo reconstruction by invoking the STAR-Fusion '--denovo_reconstruct' 
parameter. This requires that you include the '--FusionInspector inspect|validate' setting. 
Based on the STAR-Fusion predicted fusions.

FusionInspector extracts the genomic regions for the fusion partners and constructs mini-fusion-contigs containing 
the pairs of genes in their proposed fused orientation. 

fusion annotation and examination of coding effects are additionally performed on FusionInspector outputs 
if '--FusionInspector inspect|validate' is invoked in your STAR-Fusion run

--examine_coding_effect想要获得denovo组装的序列蛋白，就需要--denovo_reconstruct; --denovo_reconstruct又需要--FusionInspector inspect

输出结果：
FusionName     JunctionReadCount       SpanningFragCount       SpliceType      LeftGene        LeftBreakpoint  RightGene       RightBreakpoint LargeAnchorSupport      LeftBreakDinuc  LeftBreakEntropy        RightBreakDinuc RightBreakEntropy       FFPM
Examine Effect of Fusions on Coding Regions
Trinity-reconstructed fusion transcript
'''

class innerSoftware:
	def __init__(self):
		self.starFusion = ''
		
		

class STARFUSION:
	def __init__(self, left_fq, right_fq, genome_lib_dir, output_dir, software = innerSoftware()):
		self.left_fq = left_fq
		self.right_fq = right_fq
		self.genome_lib_dir = genome_lib_dir
		self.output_dir = output_dir
		
	def StarFusion(self):
		code = [self.software.starFusion,
			'--genome_lib_dir ' + self.genome_lib_dir,
			'--left_fq ' + self.left_fq,
			'--right_fq ' + self.right_fq,
			'--output_dir ' + self.output_dir,
			'--FusionInspector inspect',   #only the reads identified by STAR-Fusion as evidence supporting the fusion prediction are aligned directly to a target set of fusion-gene contigs for exploration using IGV
			'--denovo_reconstruct',   #If you're interested in doing de novo reconstruction of fusion transcripts based on the actual fusion-supporting RNA-Seq reads, capturing any additional variant information or novel sequence features that may be evidence among the reads, you can include de novo reconstruction by invoking the STAR-Fusion '--denovo_reconstruct' parameter.
			'--examine_coding_effect'   #explore the impact of the fusion event on coding regions
		]		
		return '   \\\n\t'.join(code)
		
		
	def printJob(self):
		return self.Star() + '\n\n'

	def __call__(self):
		subprocess.check_call(self.printJob(), shell = True)
		


def paramsParse():
	parser = argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)
	parser.add_argument('--left_fq', help = "read1 file", required = True)
	parser.add_argument('--genome_lib_dir', help = "specifies path to the genome resource lib", required = True)
	parser.add_argument('--right_fq', help = "read2 file", required = True)
	parser.add_argument('--output_dir', help = "output dir", default = '.')

	return parser.parse_args()

if __name__ == '__main__':
	args = paramsParse()
	starFusion = STARFUSION(args.left_fq, args.right_fq, args.genome_lib_dir, args.output_dir)
	starFusion()