# coding=utf-8
# guoyang
# 已完成
import subprocess
import argparse
import os.path
import sys
import os
import DepthOfCoverage
import MergeVCF

usage='''
ExomeCNV is an R package tailored to detection of CNV (Copy-Number Variants) and LOH (Loss of Heterozygosity) from exome sequencing data. It exploits the unique discrete feature of exon definitions and incredible cross-sample consistency of depth-of-coverage. ExomeCNV is most suitable when paired samples (e.g. tumor-normal pair) are available. Both of the paired samples should be processed and sequenced in a similar manner (e.g. same library prep, sequencer, average depth-of-coverage, etc.).

The inputs necessary for ExomeCNV are 1) bam/pileup file 2) exome definition file (.bed) and 3) a conservative approximation of sample admixture rate. 



The three primary steps of ExomeCNV are

    1. Calculate log coverage ratio between case and control
    2. Call CNV/LOH for each exon individually
    3. Combine exonic CNV/LOH into segments using Circular Binary Segmentation (CBS)


https://secure.genome.ucla.edu/index.php/ExomeCNV_User_Guide#Overview
###########################################

Optimization Strategy

Because there could be multiple cutoff values that satisfy min.sens and min.spec (see above), we could pick the cutoff value that optimizes for sensitivity, specificity, or an average between the two (i.e. area-under-curve, AUC = (sensitivity + specificity)/2). The choice depends on the application, but for the same reason that we should set high min.spec in the first round of calls, the user might want to optimize for "specificity" in the first round of CNV/LOH calls. In the segment-level calls, the user may elect to user "AUC" option to get the "best of both worlds," balancing between specificity and sensitivity.
It is also recommended that if the user uses the option "AUC", min.sense and min.spec should be set to the same value. By experience, when "AUC" option is used and min.sens is not equal to min.spec, the results tend to be very noisy.



Sample Admixture Rate

It is common in the case of tumor biopsy to have some heterogeneity in the sample. Particularly, some normal (non-mutated) tissues may be present in the sample at certain rate. A 30-50% admixture is not at all uncommon. The sensitivity and precision of CNV detection are dependent on a good estimate of the admixture rate, and a more conservative estimate tends to lead to fewer false positives. By default, the admixture rate "c" is set at 0.3. If the results appear to be noisy, users may consider increasing "c" up to 0.5 or more.



Circular Binary Segmentation (CBS)

In the functions multi.CNV.analyze and multi.LOH.analyze, ExomeCNV performs segmentation using CBS (implemented in DNAcopy package). There are a few parameters for CBS controlling granularity of segmentation. Intuitively, a more granulated segmentation will produce more of smaller segments. The two parameters that controls granularity of CBS are:
    sd.undo: the number of standard deviations between means to keep a breakpoint
        The higher sd.undo is, the coarser the segments will become. The default values are set to 1 and 2.
    alpha: significance levels for the test to accept change-points
        The lower alpha is, the coarse the segments will become (harder to accept a breakpoint). The default values are set to 0.05 and 0.01.


###########################################################
LOH
Choice of Test Statistics
Three test statistics based on BAF

Statistical tests that can be used in calling LOH are based on three test statistics:

    BAF as count statistic
    Variance of BAF, reflecting the amount of deviation of BAF away from its central value (~0.5)
    Absolute deviation of BAF from the null value of 0.5
    Difference between BAF's in case and control samples

Each test statistic allows for different tests and is based on different assumptions.
BAF as count statistic

These tests are only available when calling LOH on a single position (using LOH.analyze).

Options "only.tumor" and "only.normal" use only one sample (case or control) to perform binomial test against null p=0.5. We can model LOH as a binomial event, asking among N reads mapped to the position, how likely is it to observe a certain number of B-allele (BAF).

Options "two.sample.fisher" and "two.sample.prop" are similar to the binomial test for one sample above but instead of testing the observed proportion against the null value of 0.5, they compare the observed proportion between case and control. This can be modeled by binomial distribution (two.sample.prop) or hypergeometric distribution (Fisher's exact test; two.sample.fisher), hence the two possible tests.
Variance of BAF

Option "variance.f" performs F-test to compare variances of case and control BAF's
Absolute Deviation of BAF from 0.5

Options "deviation.wilcox" and "deviation.t" perform t-test and Wilcoxon Rank Sum (Mann-Whitney) Test, respectively. This is to compare the mean value of the absolute deviation of BAF from 0.5 (i.e. |BAF - 0.5|).
Difference between case and control BAF's

Option "deviation.half.norm" is based on the observation that the distribution of BAF difference between case and control are normally distributed around 0. Thus the absolute value follows folded-normal distribution. Under LOH, the absolute difference will have a higher mean value, and we can measure and test the increase in the difference using half-normal distribution.
Other Test
Cochran-Mantel-Haenszel Statistics

Option "CMH" or "mantelhaen" uses Cochran-Mantel-Haenszel Chi-sq test for common odds ratio equal to 1. It requires that the number of stata N >= 2. In case N = 1, it is equivalent to Pearson's Chi-sq (prop.test). This is useful when trying to call LOH for segments, which contain multiple heterozygous positions, each with its own contigency table. The only problem with this test is that it requires phasing information, which does not always exist. Thus it is not recommended for use. 



##################################################
The two primary steps for LOH calling are:

    Calling LOH at each heterozygous position
    Combine multiple positions into LOH segments

##################################################
Q: My data is paired end data. The ends are 50bp and 35bp. The classify.eCNV step asks for read lengths. Should I put 50, 35, or 85?

A: 42. Read length is used to convert depth-of-coverage back into approximate read counts. As long as the number of reads from both ends are roughly equal, using the mean read length should give the right approximate.

Q: I have multiple (unpaired) case and control exomes. How can I use ExomeCNV to detect CNV?

A: ExomeCNV was developed for closely matched paired samples and performs best in that setting. However, one may pool together a set of samples and compare pooled case and pooled control exomes. The function pool.coverage in ExomeCNV package is available for this use. It is important to note that the power calculation might be off and a conservative evaluation of the results is recommended.

Q: I noticed that when running function "read.coverage.gatk(file)" to reformat GATK output, these following two fields "sequenced.base" and "bease.with..10.coveratge" are all with "NA". Do you think "NA" in these two fields will affect the detection results?

A: No, the fields "sequenced.base" and "bases.with..10.coverage" are degenerate and are not important for ExomeCNV to function. They are there for a historical reason. 


应用范围：成对WES
'''

class innerSoftware:
	def __init__(self):
		self.R_2_15_3 = 'Rscript'
		self.perl = 'perl'



class EXOMECNV:
	def __init__(self, caseBam, controlBam, bed, admix, loh, caseVcf, controlVcf, outdir, software=innerSoftware()):
		self.caseBam = os.path.abspath(caseBam)
		self.controlBam = os.path.abspath(controlBam)
		self.bed = bed
		self.admix = admix
		self.loh = loh
		self.caseVcf = caseVcf
		self.controlVcf = controlVcf
		self.outdir = os.path.abspath(outdir)
		self.software = software

	def DepthOfCoverage(self):
		caseCoverage = DepthOfCoverage.DEPTHOFCOVERAGE(self.caseBam, self.bed, True, True, self.outdir)
		controlCoverage = DepthOfCoverage.DEPTHOFCOVERAGE(self.controlBam, self.bed, True, True, self.outdir)
		code = caseCoverage.printJob() + '\n' + controlCoverage.printJob()
		self.caseCoverageFile = caseCoverage.outfile + '.sample_interval_summary'
		self.controlCoverageFile = controlCoverage.outfile + '.sample_interval_summary'
		return code
		
	def ExomeCNV(self):
		self.software.exomeCNV = os.path.join(sys.path[0], 'exomeCNV.R')
		code = [self.software.R_2_15_3, self.software.exomeCNV, self.caseCoverageFile, self.controlCoverageFile, os.path.basename(self.caseBam).replace('.bam', '.somatic'), self.admix]
		return ' '.join(code)

	def CombineVariants(self):
		self.combinedVcf = os.path.join(self.outdir, os.path.basename(self.caseVcf).replace('.vcf', '.combined.vcf'))
		combineVariants = MergeVCF.MERGEVCF(self.caseVcf + ',' + self.controlVcf, self.combinedVcf, 'horizontal')
		return combineVariants.printJob()
	
	def ForLohFile(self):
		'''
		from GATK VCF file, you can use the perl script for.loh.files.pl to create BAF file
		perl for.loh.files.pl <input vcf> <output normal.baf.txt> <output tumor.baf.txt> <normal_col_#> <tumor_col_#>
		'''
		self.software.forlohfile = os.path.join(sys.path[0], 'for.loh.files.pl')
		self.caseBaf = os.path.join(self.outdir, os.path.basename(self.caseVcf).replace('.vcf', '.baf'))
		self.controlBaf = os.path.join(self.outdir, os.path.basename(self.controlVcf).replace('.vcf', '.baf'))
		code = [self.software.perl, self.softwore.forlohfile, self.combinedVcf, self.controlBaf, self.caseBaf, '2', '1']
		return ' '.join(code)

	def ExomeLOH(self):
		self.software.exomeLOH = os.path.join(sys.path[0], 'exomeLOH.R')
		code = [self.software.R_2_15_3, self.software.exomeLOH, self.caseBaf, self.controlBaf, os.path.basename(self.caseBam).replace('.bam', '.LOH')]
		return ' '.join(code)
		
	def printJob(self):
		code = [self.DepthOfCoverage(), self.ExomeCNV()]
		if self.loh:
			code += [self.CombineVariants(), self.ForLohFile(), self.ExomeLOH()]
		return ' &&\n'.join(code)

	def __call__(self):
		subprocess.check_call(self.printJob(), shell = True)

	


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--caseBam', help="case bam file",required=True)
	parser.add_argument('--controlBam', help="control bam file", required=True)
	parser.add_argument('--bed', help='exome definition file (.bed)', required=True)
	parser.add_argument('--admix', help="a conservative approximation of sample admixture rate", default=0.5, type = float)
	parser.add_argument('--loh', help="LOH calling", action='store_true')
	parser.add_argument('--caseVcf', help="snv calling vcf from case GATK", default=None)
	parser.add_argument('--controlVcf', help="snv calling vcf from control GATK", default=None)
	parser.add_argument('--outdir', help="output dir", default=".")
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	exomeCNV = EXOMECNV(args.caseBam, args.controlBam, args.bed, args.admix, args.loh, args.caseVcf, args.controlVcf, args.outdir)
	exomeCNV()