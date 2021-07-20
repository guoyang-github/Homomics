# coding=utf-8
# guoyang
#20170627

#1.基本的实现原理   2.适应的输入数据和应用情况是什么   3.能得到怎样的结果

# FusionMap can detect fusion events in both single- and paired-end datasets from either RNA-Seq or gDNA-Seq studies and characterize fusion junctions at base-pair resolution.
# FusionMap includes reference indexing, read filtering, fusion alignment and reporting in one package. 
# When the read length is long or the inner distance (the gap size) between read pairs is short, the number of junction-spanning reads would be larger than that of discordant read pairs. Then, fusion detection based on junction-spanning reads is more powerful than the paired-end approach.
# http://www.arrayserver.com/wiki/index.php?title=FusionMap
# DNA RNA 单端 双端都可以，主要是利用junction-spanning reads

# 所需的参考基因组和基因模型可以自动下载预编译好的，也可以手动下载放置在预设的路径下，也可以自行制作，方法如下：
#Building a Reference Library 
#   FusionMap.exe --buildref FusionMap_Base_Dir fasta_file_name ref_lib_name
#Building a Gene Model 
#   FusionMap.exe --buildgm FusionMap_Base_Dir gtf_file_name ref_lib_name gene_model_name

# 预编译好的reference和gene model在这里查看
# http://www.arrayserver.com/wiki/index.php?title=A_list_of_compiled_genome_and_gene_model_from_OmicSoft

#Fusion 结果过滤
#FusionMap reports as many fusion candidates as possible and provides multiple features for users to filter out false positives. To reduce false positives for REAL data, the recommended filtering sets are
#    sample.SeedCount >= 3
#    SplicePatternClass=CanonicalPattern[Major] or CanonicalPattern[Minor]
#    Filter=Empty 
#To be more stringent, user can further filter using
#    FrameShiftClass=InFrame
#    OnExonBoundary=Both 


#Fusion结果注释
#http://www.arrayserver.com/wiki/index.php?title=Fusion_SE_report


#另FusionMap也可以利用成对reads，PE关系做Fusion detection
#见Fusion PE Utility



import subprocess
import argparse
import os.path
import sys
import os

usage = '''
FusionMap.exe --semap /pathTo/FusionMap_Base_Dir Human.B37.3 RefGene /test/secontrol.txt > /test/run.log
'''

controlFile = '''
<Files>
{fq1}
{fq2}

<Options>
//MonoPath option is required when path to mono are not in PATH and job cannot start for spawn off jobs
MonoPath={mono}
PairedEnd=True			             //Automatically pair two fastq files as one sample to run fusion analysis
RnaMode={type}			             //Detect fusion results. If True, specifies that the samples are RNA-Seq. IF set to False, specifies that the samples are DNA-Seq
ThreadNumber=8			             //Possible values: 1-100. Default value]]=1
FileFormat=FASTQ		             //Possible values: FASTQ, QSEQ, FASTA. Default value]]=FASTQ
CompressionMethod=Gzip		         //Gzip formatted input files
Gzip=True			                 //Gzip
QualityEncoding=Automatic	         //Auto detect quality coding in the fastq file or specify with Illumina or Sanger
//AutoPenalty=True		             //Set alignment penalty cutoff to automatic based on read length: Max (2,(read length-31)/15)
//FixedPenalty=2			             //If AutoPenalty=False, Fixed Penalty will be used
//FilterUnlikelyFusionReads=True       //Enable filtering step
//FullLengthPenaltyProportion=8	     //Filtering normal reads allowing 8% of alignment mismatches of the reads
//MinimalFusionAlignmentLength=0	     //Default (alpha in the paper) value=0 and the program will automatically set minimal Seed Read end length to Min(20, Max(17,floor(ReadLength/3))). The program will use the specified value if user sets any > 0.
//FusionReportCutoff=1		         //# of allowed multiple hits of read ends; Possible values: 1-5. Default value=1 (beta in paper); 
//NonCanonicalSpliceJunctionPenalty=4  //Possible values: 0-10. Default value= 2 (G); 
//MinimalHit=4			             //Minimal distinct fusion read; Possible values: 1-10000, Default value=2 
//MinimalRescuedReadNumber=1	         //Minimal rescued read number. Default value= 1
//MinimalFusionSpan=5000		         //Minimal distance (bp) between two fusion breakpoints
//RealignToGenome=True		         //If True, seed read ends are re-aligned to genome to see if it is <= FusionReportCutoff in RNA-Seq.
OutputFusionReads=True		         //Out put Fusion reads as BAM files for genome browser. Default value=True

<Output>
TempPath={TempPath}
OutputPath={OutputPath}
OutputName=FusionMapReport
'''




class innerSoftware:
	def __init__(self):
		self.mono = '/mono-2.10.9/bin/mono'
		self.fusionmap = '/FusionMap_2015-03-31/bin/FusionMap.exe'
		

class FUSIONMAP:
    def __init__(self, fastq1, fastq2, outdir, type, ref, geneModel, software=innerSoftware()):
        self.fastq1 = os.path.abspath(fastq1)
        self.fastq2 = os.path.abspath(fastq2)
        self.outdir = os.path.abspath(outdir)
        self.type = type
        self.ref = ref
        self.geneModel = geneModel
        self.software=software
		
    def FusionMap(self):
        self.controlFile = os.path.join(self.outdir, 'fusionmap.controlfile')
        self.TempPath = os.path.join(self.outdir, 'TMP')
        with open(self.controlFile, 'w') as controlF:
            controlF.write(controlFile.format(mono = self.software.mono, fq1 = self.fastq1, fq2 = self.fastq2, type = self.type, TempPath = self.TempPath, OutputPath = self.outdir))
        code = [self.software.mono,
                self.software.fusionmap,
                '--semap',
                self.outdir,
                self.ref,
                self.geneModel,
                self.controlFile,
                '>',
                os.path.join(self.outdir, 'run.log')]
        return ' '.join(code)
		
		
    def printJob(self):
        return self.FusionMap() + '\n\n'

    def __call__(self):
        subprocess.check_call(self.FusionMap(), shell=True)
		


def paramsParse():
    parser=argparse.ArgumentParser(description = usage, formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument('--fastq1', help = "read 1 fastq", required = True)
    parser.add_argument('--fastq2', help = "read 2 fastq", required = True)
    parser.add_argument('--outdir', help = "output dir", default = '.')
    parser.add_argument('--type', help = 'RNA or DNA', default = 'DNA', choices = ['RNA', 'DNA'])
    parser.add_argument('--ref', help = 'a compiled genome.', default = 'Human.B37.3')
    parser.add_argument('--geneModel', help = 'gene model', default = 'RefGene')
    return parser.parse_args()

if __name__ == '__main__':
    args=paramsParse()
    fusionMap=FUSIONMAP(args.fastq1, args.fastq2, args.outdir, args.type, args.ref, args.geneModel)
    fusionMap()