# coding=utf-8
# guoyang

import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar picard.jar CollectWgsMetrics \
       I=input.bam \
       O=collect_wgs_metrics.txt \
       R=reference_sequence.fasta 


Metrics for evaluating the performance of whole genome sequencing experiments.
Field	            Description
GENOME_TERRITORY	The number of non-N bases in the genome reference over which coverage will be evaluated.
MEAN_COVERAGE	    The mean coverage in bases of the genome territory, after all filters are applied.
SD_COVERAGE	        The standard deviation of coverage of the genome after all filters are applied.
MEDIAN_COVERAGE	    The median coverage in bases of the genome territory, after all filters are applied.
MAD_COVERAGE	    The median absolute deviation of coverage of the genome after all filters are applied.
PCT_EXC_MAPQ	    The fraction of aligned bases that were filtered out because they were in reads with low mapping quality (default is < 20).
PCT_EXC_DUPE	    The fraction of aligned bases that were filtered out because they were in reads marked as duplicates.
PCT_EXC_UNPAIRED	The fraction of aligned bases that were filtered out because they were in reads without a mapped mate pair.
PCT_EXC_BASEQ	    The fraction of aligned bases that were filtered out because they were of low base quality (default is < 20).
PCT_EXC_OVERLAP	    The fraction of aligned bases that were filtered out because they were the second observation from an insert with overlapping reads.
PCT_EXC_CAPPED	    The fraction of aligned bases that were filtered out because they would have raised coverage above the capped value (default cap = 250x).
PCT_EXC_TOTAL	    The total fraction of aligned bases excluded due to all filters.
PCT_1X	            The fraction of bases that attained at least 1X sequence coverage in post-filtering bases.
PCT_5X	            The fraction of bases that attained at least 5X sequence coverage in post-filtering bases.
PCT_10X	            The fraction of bases that attained at least 10X sequence coverage in post-filtering bases.
PCT_15X	            The fraction of bases that attained at least 15X sequence coverage in post-filtering bases.
PCT_20X	            The fraction of bases that attained at least 20X sequence coverage in post-filtering bases.
PCT_25X	            The fraction of bases that attained at least 25X sequence coverage in post-filtering bases.
PCT_30X	            The fraction of bases that attained at least 30X sequence coverage in post-filtering bases.
PCT_40X	            The fraction of bases that attained at least 40X sequence coverage in post-filtering bases.
PCT_50X	            The fraction of bases that attained at least 50X sequence coverage in post-filtering bases.
PCT_60X	            The fraction of bases that attained at least 60X sequence coverage in post-filtering bases.
PCT_70X	            The fraction of bases that attained at least 70X sequence coverage in post-filtering bases.
PCT_80X	            The fraction of bases that attained at least 80X sequence coverage in post-filtering bases.
PCT_90X	            The fraction of bases that attained at least 90X sequence coverage in post-filtering bases.
PCT_100X	        The fraction of bases that attained at least 100X sequence coverage in post-filtering bases.
HET_SNP_SENSITIVITY	The theoretical HET SNP sensitivity.
HET_SNP_Q	        The Phred Scaled Q Score of the theoretical HET SNP sensitivity.

'''

class innerSoftware:
	def __init__(self):
		self.java1_8='java'
		self.picard='picard.jar'

class innerDatabase:
	def __init__(self):
		self.reference='human_g1k_v37_decoy.fasta'
		self.referenceIndex='human_g1k_v37_decoy.fasta.fai'

class WGSMETRICS:
	def __init__(self, infile, outdir, interval, software=innerSoftware(), database=innerDatabase()):
		self.outdir=os.path.abspath(outdir)
		self.infile=os.path.abspath(infile)
		self.interval=interval
		self.software=software
		self.database=database
		
	def WgsMetrics(self):
		self.outfile=os.path.join(self.outdir, os.path.basename(self.infile) + '_collect_wgs_metrics.txt')
		code=[self.software.java1_8 + ' -jar ' + self.software.picard + ' CollectWgsMetrics',
			'INPUT=' + self.infile,
			'OUTPUT=' + self.outfile,
			'REFERENCE_SEQUENCE=' + self.database.reference,
			#'MINIMUM_MAPPING_QUALITY='   #Minimum mapping quality for a read to contribute coverage. Default value: 20. This option can be set to 'null' to clear the default value. 
			#'MINIMUM_BASE_QUALITY='   #Minimum base quality for a base to contribute coverage. N bases will be treated as having a base quality of negative infinity and will therefore be excluded from coverage regardless of the value of this parameter. Default value: 20. This option can be set to 'null' to clear the default value.
			#'COVERAGE_CAP='   #Treat positions with coverage exceeding this value as if they had coverage at this value (but calculate the difference for PCT_EXC_CAPPED). Default value: 250. This option can be set to 'null' to clear the default value. 
			'INCLUDE_BQ_HISTOGRAM=true',   #Determines whether to include the base quality histogram in the metrics file. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false} 
			#'COUNT_UNPAIRED=',   #If true, count unpaired reads, and paired reads with one end unmapped Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false} 
			#'SAMPLE_SIZE='    #Sample Size used for Theoretical Het Sensitivity sampling. Default is 10000. Default value: 10000. This option can be set to 'null' to clear the default value. 
		]		
		if self.interval:
			code+=['INTERVALS=' + self.interval]   #An interval list file that contains the positions to restrict the assessment. Please note that all bases of reads that overlap these intervals will be considered, even if some of those bases extend beyond the boundaries of the interval. The ideal use case for this argument is to use it to restrict the calculation to a subset of (whole) contigs. To restrict the calculation just to individual positions without overlap, please see CollectWgsMetricsFromSampledSites. Default value: null. 
		return '   \\\n\t'.join(code)
		
	def __call__(self):
		return self.WgsMetrics() + '\n\n'

	def runJob(self):
		subprocess.check_call(self.WgsMetrics(), shell=True)
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--infile', help="bam file",required=True)
	parser.add_argument('--outdir', help="output dir", default='.')
	parser.add_argument('--interval', help="An interval list file that contains the positions to restrict the assessment. ", default=None)
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	wgsMetrics=WGSMETRICS(args.infile, args.outdir, args.interval)
	wgsMetrics()