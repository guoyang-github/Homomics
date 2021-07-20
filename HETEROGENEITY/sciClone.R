#Rscript sciClone.R vaf exclude cn sample outdir
#sciClone integrates the read depth and copy number information at single nucleotide variant locations and clusters the variants in copy neutral regions, to formalize description of the sub-clonal architecture of the sample.

#sciClone(vafs, copyNumberCalls=NULL, regionsToExclude=NULL,
#                  sampleNames, minimumDepth=100, clusterMethod="bmm",
#                  clusterParams=NULL, cnCallsAreLog2=FALSE,
#                  useSexChrs=TRUE, doClustering=TRUE, verbose=TRUE,
#                  copyNumberMargins=0.25, maximumClusters=10,
#                  annotation=NULL, doClusteringAlongMargins=TRUE,
#                  plotIntermediateResults=0)

Args <- commandArgs(trailingOnly = TRUE)
#获得参数
vaf <- strsplit(Args[1], split=',')
volume <- length(vaf)
exclude <- Args[2] #or 'NULL'
cn <- strsplit(Args[3], split=',')
if (length(cn) != volume) stop('file number error')
sample <- strsplit(Args[4], split=',')
if (length(sample) != volume) stop('file number error')
outdir <- Args[5]


library(sciClone)
#整理参数
vafs <- list()
cns <- list()
for (i in 1:volume)
	vafs[i] <- read.table(vaf[i], header=T)
	cns[i] <- read.table(cn[i])
if (exclude != 'NULL')	
	regions <- read.table(exclude)
esle
	regions <- NULL	
names <- sample

#程序主体
sc <- sciClone(vafs=vafs, copyNumberCalls=cns, sampleNames=names, regionsToExclude=regions，minimumDepth=8)
#vafs: a list of dataframes containing variant allele fraction data for single nucleotide variants in 5-column format: 1.chromosome 2. position 3. reference-supporting read counts 4. variant-supporting read counts 5. variant allele fraction(between 0-100)
#copyNumberCalls: list of dataframes containing copy number segments in 4-column format: 1. chromosome 2. segment start position 3.segment stop position 4. copy number value for that segment. Unrepresented regions are assumed to have a copy number of 2.
#regionsToExclude: Exclusion regions in 3-column format: 1. chromosome 2. window start position 3. window stop position; Single nucleotide variants falling into these windows will not be included in the analysis. Use this input for LOH regions, for example.
#sampleNames: vector of names describing each sample ex: ("Primary", Tumor", "Relapse")
#minimumDepth: threshold used for excluding low-depth variants
#maximumClusters: max number of clusters to consider when choosing the component fit to the data.
#cnCallsAreLog2: boolean argument specifying whether or not the copy number predictions are in log2 format (as opposed to being absolute copy number designations)
#copyNumberMargins: In order to identify cleanly copy-number neutral regions, sciClone only considers sites with a copy number of 2.0 +/- this value. For example, if set to 0.25, regions at 2.20 will be considered cn-neutral, and regions at, 2.30 will not.

writeClusterTable(sc, paste(outdir, '/', cluster, volume))
sc.plot1d(sc, paste(outdir, '/', cluster, volume, '.1d.pdf'))
if (volume >= 2) sc.plot2d(sc, paste(outdir, '/', cluster, volume, '.2d.pdf'))
if (volume == 3) sc.plot3d(sc, sc@sampleNames, size=700, outputFile=paste(outdir, '/', cluster, volume, '.3d.gif'))
