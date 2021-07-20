#https://secure.genome.ucla.edu/index.php/ExomeCNV_User_Guide#Overview  ExomeCNV_User_Guide

#The three primary steps of ExomeCNV are

#    1. Calculate log coverage ratio between case and control
#    2. Call CNV/LOH for each exon individually
#    3. Combine exonic CNV/LOH into segments using Circular Binary Segmentation (CBS)


#The two primary steps for LOH calling are:

#    Calling LOH at each heterozygous position
#    Combine multiple positions into LOH segments


library(ExomeCNV)
library(DNAcopy)


Args <- commandArgs(trailingOnly = TRUE)
tumorBaf = Args[1]
normalBaf = Args[2]
name = Args[3]


#Calling LOH on each heterozygous position
eLOH = LOH.analyze(normalBaf, tumorBaf, alpha=0.05, method="two.sample.fisher")

#Combine multiple positions into LOH segments
the.loh = multi.LOH.analyze(normalBaf, tumorBaf, all.loh.ls=list(eLOH), test.alpha=0.001, method="variance.f", sdundo=c(0,0), alpha=c(0.05,0.01))
do.plot.loh(the.loh, normalBaf, tumorBaf, "two.sample.fisher", plot.style="baf")
write.loh.output(the.loh, paste(name,".eloh",sep=''))

#This will give LOH intervals in a format similar to a bed file. If you wish to assign LOH status to all heterozygous positions, expand.loh function may be used: 
expanded.loh = expand.loh(the.loh, normalBaf)
write.loh.output(expanded.loh, paste(name, ".all", sep=''))
