# coding=utf-8
# guoyang
# 封装软件  1. 独立执行  2.插入流程
import subprocess
import argparse
import os.path
import sys
import os

usage='''
 java -jar GenomeAnalysisTK.jar \
   -T DepthOfCoverage \
   -R reference.fasta \
   -o file_name_base \
   -I input_bams.list         这里是一个bam文件的list，可以是多个文件
   [-geneList refSeq.sorted.txt] \
   [-pt readgroup] \
   [-ct 4 -ct 6 -ct 10] \
   [-L my_capture_genes.interval_list]

This tool processes a set of bam files to determine coverage at different levels of partitioning and aggregation. Coverage can be analyzed per locus, per interval, per gene, or in total; can be partitioned by sample, by read group, by technology, by center, or by library; and can be summarized by mean, median, quartiles, and/or percentage of bases covered to or beyond a threshold. Additionally, reads and bases can be filtered by mapping or base quality score. 



OUTPUT:
 Tables pertaining to different coverage summaries. Suffix on the table files declares the contents:

    no suffix: per locus coverage
    _summary: total, mean, median, quartiles, and threshold proportions, aggregated over all bases
    _statistics: coverage histograms (# locus with X coverage), aggregated over all bases
    _interval_summary: total, mean, median, quartiles, and threshold proportions, aggregated per interval
    _interval_statistics: 2x2 table of # of intervals covered to >= X depth in >=Y samples
    _gene_summary: total, mean, median, quartiles, and threshold proportions, aggregated per gene
    _gene_statistics: 2x2 table of # of genes covered to >= X depth in >= Y samples
    _cumulative_coverage_counts: coverage histograms (# locus with >= X coverage), aggregated over all bases
    _cumulative_coverage_proportions: proprotions of loci with >= X coverage, aggregated over all bases


关于dup的处理：
Reads that are marked in the bam as duplicates are removed by the Depth Of Coverage tool automatically - there's nothing you need to do, assuming you have used a tool to mark the duplicates first (see our best practices for more details).
http://gatkforums.broadinstitute.org/gatk/discussion/1392/remove-duplicate-reads-from-the-depthofcoverage

'''

class innerSoftware:
    def __init__(self):
        self.java1_8='java'
        self.gatk='GenomeAnalysisTK.jar'

class innerDatabase:
    def __init__(self):
        self.reference='human_g1k_v37_decoy.fasta'
        self.referenceIndex='human_g1k_v37_decoy.fasta.fai'


class DEPTHOFCOVERAGE:
    def __init__(self, infile, interval, omitBaseOutput, omitLocusTable, outdir, software=innerSoftware(), database=innerDatabase()):
        self.infile=os.path.abspath(infile)   #The GATK reads argument (-I, --input_file) supports only BAM/CRAM files with the .bam/.cram extension and lists of BAM/CRAM files with the .list extension
        self.interval=interval
        self.omitBaseOutput=omitBaseOutput
        self.omitLocusTable=omitLocusTable
        self.outdir=os.path.abspath(outdir)
        self.software=software
        self.database=database
		
    def DepthOfCoverage(self):
        self.outfile=os.path.join(self.outdir, os.path.basename(self.infile) + '.DOC')
        code=[self.software.java1_8 + ' -jar ' + self.software.gatk,
            '-T DepthOfCoverage',
            '-R ' + self.database.reference,
            '--out ' + self.outfile,   #file_name_base
            '-I ' + self.infile,
			#'--calculateCoverageOverGenes',   #Calculate coverage statistics over this list of genes  Specify a RefSeq file for use in aggregating coverage statistics over genes. This argument is incompatible with --calculateCoverageOverGenes and --omitIntervalStatistics. A warning will be logged and no output file will be produced for the gene list if these arguments are enabled together.
            #'--countType COUNT_READS',   #How should overlapping reads from the same fragment be handled?
                #COUNT_READS
                #    Count all reads independently (even if from the same fragment).
                #COUNT_FRAGMENTS
                #    Count all fragments (even if the reads that compose the fragment are not consistent at that base).
                #COUNT_FRAGMENTS_REQUIRE_SAME_BASE
                #    Count all fragments (but only if the reads that compose the fragment are consistent at that base). 
            #'--maxBaseQuality',
            #'--maxMappingQuality',
            #'--minBaseQuality',
            #'--minMappingQuality',
            '--outputFormat rtable',   #Output file format (e.g. csv, table, rtable); defaults to r-readable table.
            '--partitionType sample',   #Partition type for depth of coverage  By default, coverage is partitioning by sample, but it can be any combination of sample, readgroup and/or library
            '--summaryCoverageThreshold 1',   #Coverage threshold (in percent) for summarizing statistics  For summary file outputs, report the percentage of bases covered to an amount equal to or greater than this number (e.g. % bases >= CT for each sample). Defaults to 15; can take multiple arguments.
            '--summaryCoverageThreshold 4',
            '--summaryCoverageThreshold 10',

            #'--omitDepthOutputAtEachBase',   #Do not output depth of coverage at each base   --omitBaseOutput
            #'--omitIntervalStatistics false',   #Do not calculate per-interval statistics
            #'--omitLocusTable',   #Do not calculate per-sample per-depth counts of loci
            #'--omitPerSampleStats',   #Do not output the summary files per-sample
            #'--printBaseCounts',   #Add base counts to per-locus output  Instead of reporting depth, the program will report the base pileup at each locus 
            #'--ignoreDeletionSites',
            #'--includeDeletions',
            #'--includeRefNSites'			
		]
        if self.omitBaseOutput:
            code+=['--omitDepthOutputAtEachBase']   #Disabling the tabulation of total coverage at every base should speed up processing.
        if self.omitLocusTable:
            code+=['--omitLocusTable']   #Disabling the tabulation of locus statistics (# loci covered by sample by coverage) should speed up processing. 
        if self.interval:
            code+=['-L ' + self.interval]
        return '   \\\n\t'.join(code)
		
		
    def printJob(self):
        return self.DepthOfCoverage() + '\n\n'

    def __call__(self):
        subprocess.check_call(self.DepthOfCoverage(), shell=True)
		


def paramsParse():
    parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--infile',help="bam file list with bai index in the same dir",required=True)
    parser.add_argument('--interval', help="capture regions", default=None)
    parser.add_argument('--omitBaseOutput', help="Do not output depth of coverage at each base", action='store_true')
    parser.add_argument('--omitLocusTable', help="Do not calculate per-sample per-depth counts of loci", action='store_true')
    parser.add_argument('--outdir', help="output dir", default='.')
    return parser.parse_args()

if __name__ == '__main__':
    args=paramsParse()
    depthOfCoverage=DEPTHOFCOVERAGE(args.infile, args.interval, args.omitBaseOutput, args.omitLocusTable, args.outdir)
    depthOfCoverage()