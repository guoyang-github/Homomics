#https://secure.genome.ucla.edu/index.php/ExomeCNV_User_Guide#Overview  ExomeCNV_User_Guide

#The three primary steps of ExomeCNV are

#    1. Calculate log coverage ratio between case and control
#    2. Call CNV/LOH for each exon individually
#    3. Combine exonic CNV/LOH into segments using Circular Binary Segmentation (CBS)


Args <- commandArgs(trailingOnly = TRUE)
tumorCoverage = Args[1]
normalCoverage = Args[2]
name = Args[3]
admix = Args[4]


chr.list.ucsc.hg19 = c(paste('chr',1:22,sep=''), 'chrX', 'chrY')  #default
chr.list.1000g.b37 = c(as.character(1:22), 'X', 'Y')
chr.list = chr.list.ucsc.hg19



#Calculate log coverage ratio
library(ExomeCNV)
tumor = read.coverage.gatk(tumorCoverage)
normal = read.coverage.gatk(normalCoverage)
coverageRatio.logR = calculate.logR(normal, tumor)


#Call CNV for each exon
#Then call CNV on each exon (using classify.eCNV), one chromosome at a time. We recommend high min.spec (0.9999) and option="spec" to be conservative against false positive. This is because whatever is called at exon level will persist through merging step (Step 3). 
aggregation.eCNV = c()
for (i in 1:length(chr.list)) {
	idx = (normal$chr == chr.list[i])
	ecnv = classify.eCNV(normal=normal[idx,], tumor=tumor[idx,], logR=coverageRatio.logR[idx], min.spec=0.9999, min.sens=0.9999, option="spec", admix=admix, read.len=150)
	#classify.eCNV: Calculate specificity and sensitivity (power) of detecting CNV based on depth of coverage and log ratio of all exons. Make a call when sufficient specificity and sensitivity are achieved.
	#admix: contamination rate (admixture rate), the proportion of the normal cells in the tumor samples.
	#read.len: sequence read length.
	#option: objective quantity to optimize over when minimum sensitivity and specificity are achieved. Possible opetions are ‘sens’ for sensitivity, ‘spec’ for specificity, ‘auc’ for area under curve = (specificity + sensitivity)/2.
	#test.num.copy： copy numbers to be tested. 1 for deletion, 3 for duplication, 4 and beyond for amplification. Default to (1,3,4,5)
	aggregation.eCNV = rbind(aggregation.eCNV, ecnv)
}


#Combine exonic CNV into larger segments
#Here, we use lower min.spec and min.sens and option="auc" to be less conservative and allow for more discovery. 
library(DNAcopy)
results.cnv = multi.CNV.analyze(normal, tumor, logR=coverageRatio.logR, all.cnv.ls=list(aggregation.eCNV), coverage.cutoff=5, min.spec=0.99, min.sens=0.99, option="auc", admix=admix)


#From here we can plot the results and export outputs
do.plot.eCNV(aggregation.eCNV, lim.quantile=0.99, style="idx", line.plot=F)
do.plot.eCNV(results.cnv, lim.quantile=0.99, style="bp", bg.cnv=aggregation.eCNV, line.plot=T)
#Generate output files from ExomeCNV outputs. The files geneated are: 1. .cnv.txt file with all CNV calls 2. .exon.lrr.txt file containing log coverage ratio for each exon 3. .segment.lrr.txt file containing log coverage ratio for each segment (as defined by CBS) 4. .segment.copynumber.txt file containing copy number calls for each segment 5. .cnv.png file, a plot of the results
write.output(aggregation.eCNV, results.cnv, name)  #需要在源代码包里面给这个函数的png函数增加cairo库，再安装

