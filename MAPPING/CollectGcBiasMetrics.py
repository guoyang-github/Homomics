# coding=utf-8
# guoyang
import subprocess
import argparse
import os.path
import sys
import os

usage='''
java -jar picard.jar CollectGcBiasMetrics \
      I=input.bam \
      O=gc_bias_metrics.txt \
      CHART=gc_bias_metrics.pdf \
      S=summary_metrics.txt \
      R=reference_sequence.fasta

Detailed metrics：
(GC）The table of detailed metrics includes GC percentages for each bin
(WINDOWS) the percentage of WINDOWS corresponding to each GC bin of the reference sequence.
The mean of the distribution will vary among organisms; human DNA has a mean GC content of 40%, suggesting a slight preponderance of AT-rich regions.
(READ_STARTS) the numbers of reads that start within a particular %GC content bin (READ_STARTS)
(MEAN_BASE_QUALITY) and the mean base quality of the reads that correspond to a specific GC content distribution window (MEAN_BASE_QUALITY). 
(NORMALIZED_COVERAGE) NORMALIZED_COVERAGE is a relative measure of sequence coverage by the reads at a particular GC content.


SUMMARY TABLE:
Field	Description
ACCUMULATION_LEVEL	
READS_USED  This option is used to mark including or excluding duplicates.
WINDOW_SIZE The window size on the genome used to calculate the GC of the sequence.
TOTAL_CLUSTERS	The total number of clusters that were seen in the gc bias calculation.
ALIGNED_READS	The total number of aligned reads used to compute the gc bias metrics.
AT_DROPOUT	Illumina-style AT dropout metric. Calculated by taking each GC bin independently and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[0..50].
GC_DROPOUT	Illumina-style GC dropout metric. Calculated by taking each GC bin independently and calculating (%ref_at_gc - %reads_at_gc) and summing all positive values for GC=[50..100].
which indicate the percentage of misaligned reads that correlate with low (%-GC is < 50%) or high (%-GC is > 50%) GC content respectively.

GC_NC_0_19	Normalized coverage over quintile of GC content ranging from 0 - 19.
GC_NC_20_39	Normalized coverage over each quintile of GC content ranging from 20 - 39.
GC_NC_40_59	Normalized coverage over each quintile of GC content ranging from 40 - 59.
GC_NC_60_79	Normalized coverage over each quintile of GC content ranging from 60 - 79.
GC_NC_80_100	Normalized coverage over each quintile of GC content ranging from 80 - 100.
The percentage of 'coverage' or depth in a GC bin is calculated by dividing the number of reads of a particular GC content by the mean number of reads of all GC bins. A number of 1 represents mean coverage, a number less than 1 represents lower than mean coverage (e.g. 0.5 means half as much coverage as average) while a number greater than 1 represents higher than mean coverage (e.g. 3.1 means this GC bin has 3.1 times more reads per window than average). This tool also tracks mean base-quality scores of the reads within each GC content bin, enabling the user to determine how base quality scores vary with GC content. 


图像：
The chart output associated with this data table plots the NORMALIZED_COVERAGE, the distribution of WINDOWs corresponding to GC percentages, and base qualities corresponding to each %GC bin.
红色的就reference GC分布的本质属性
绿色是质量值  This tool also tracks mean base-quality scores of the reads within each GC content bin, enabling the user to determine how base quality scores vary with GC content.
蓝色coverage随reference上GC含量的分布




针对WGS，用做WES应该是不合适的
'''

class innerSoftware:
    def __init__(self):
        self.java1_8='java'
        self.picard='picard.jar'

class innerDatabase:
    def __init__(self):
        self.reference='human_g1k_v37_decoy.fasta'

class GCBIASMETRICS:
    def __init__(self, infile, outdir, software=innerSoftware(), database=innerDatabase()):
        self.outdir=os.path.abspath(outdir)
        self.infile=os.path.abspath(infile)
        self.software=software
        self.database=database
		
    def GcBiasMetrics(self):
        self.outfile=os.path.join(self.outdir, os.path.basename(self.infile) + '_gc_bias_metrics.txt')
        code=[self.software.java1_8 + ' -jar ' + self.software.picard + ' CollectGcBiasMetrics',
            'INPUT=' + self.infile,
            'OUTPUT=' + self.outfile,
            'CHART_OUTPUT=' + self.outfile + '.pdf',   #The PDF file to render the chart to. Required. 
            'SUMMARY_OUTPUT=' + self.outfile + '_summary.txt',   #The text file to write summary metrics to. Required. 
            'REFERENCE_SEQUENCE=' + self.database.reference,
            #'ALSO_IGNORE_DUPLICATES=true'   #这个参数新版picard好像没有了Use to get additional results without duplicates. This option allows to gain two plots per level at the same time: one is the usual one and the other excludes duplicates. Default value: false.
            #'MINIMUM_GENOME_FRACTION'   #For summary metrics, exclude GC windows that include less than this fraction of the genome. Default value: 1.0E-5. 
            #'SCAN_WINDOW_SIZE'   #The size of the scanning windows on the reference genome that are used to bin reads. Default value: 100. 
            #'IS_BISULFITE_SEQUENCED'   #Whether the SAM or BAM file consists of bisulfite sequenced reads. Default value: false. 
            #'METRIC_ACCUMULATION_LEVEL'   #The level(s) at which to accumulate metrics.  Default value: [ALL_READS]. This option can be set to 'null' to clear the default value. Possible values: {ALL_READS, SAMPLE, LIBRARY, READ_GROUP} This option may be specified 0 or more times. This option can be set to 'null' to clear the default list.
            #'ASSUME_SORTED'   #If true (default), then the sort order in the header file will be ignored.  Default value: true.
            ]		
        return '   \\\n\t'.join(code)
		
    def __call__(self):
        return self.GcBiasMetrics() + '\n\n'

    def runJob(self):
        subprocess.check_call(self.GcBiasMetrics(), shell=True)
		


def paramsParse():
    parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--infile', help="Input SAM or BAM file sorted. Required.",required=True)
    parser.add_argument('--outdir', help="DIR to write the output to", default='.')
    return parser.parse_args()

if __name__ == '__main__':
    args=paramsParse()
    gcbiasMetrics=GCBIASMETRICS(args.infile, args.outdir)
    gcbiasMetrics.runJob()