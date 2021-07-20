library(optparse)

option.list <- list(
    make_option("--snv", help = "transfered from varscan2"),
    make_option("--cnv", help = "transfered from controlfreec"),
    make_option("--region", help = "Regional boundary for mutations included during clustering.")
    make_option("--outdir", help = "output dir", default = '.')
    )
#Parameter region:
#Matrix in which each row corresponds to a genomic segment.
#Columns must include:
    #*chr* - the chromosome of the segment;
    #*start* - the first genomic position of the segment;
    #*end* - the last genomic position of the segment.
    #Default: SureSelectExome_hg19, comprising ca. 468 MB centered
    #on the human exome. Alternative user supplied regions should
    #also be coding regions, as the seletive pressure is higher as
    #compared to non-coding regions.
opt <- parse_args(OptionParser(option_list=option.list))


setwd(opt[['outdir']])
snv <- as.matrix(read.table(opt[['snv']], header = T))
cnv <- as.matrix(read.table(opt[['cnv']], header = T))
region <- as.matrix(read.table(opt[['region']], header = F))

#parameters
max_PM = 6   #Upper threshold for the number of amplicons per mutated cell. Increasing the value of this variable is not recommended unless extensive depth and breadth of coverage underly the measurements of copy numbers and allele frequencies. See also 'cellfrequency_pdf.'
min_CellFreq = 0.1   #Lower boundary for the cellular prevalence interval of a mutated cell. Mutations for which allele frequency * copy number are below min_{CellFreq}, are excluded from further computation. Decreasing the value of this variable is not recommended unless extensive depth and breadth of coverage underly the measurements of copy numbers and allele frequencies.
maxScore = 2   #Upper threshold for the noise score of subpopulation detection. Only subpopulations identified at a score below maxScore are kept.
precision = 0.018   #Precision with which subpopulation size is predicted, a small value reflects a high resolution and can lead to a higher number of predicted subpopulations.


out = runExPANdS(snv, cbs, region = region)

#结果List with fields:
#finalSPs: Matrix of predicted subpopulations. Each row corresponds to a
#          subpopulation and each column contains information about that
#          subpopulation, such as the size in the sequenced tumor bulk
#          (column *Mean Weighted*) and the noise score at which the
#          subpopulation has been detected (column *score*).
#
#      dm: Matrix containing the input mutations with at least seven
#          additional columns:
#          *SP* - the subpopulation to which the point mutation has been
#          asssigned;
#          *SP_cnv* - the subpopulation to which the CNV has been
#          asssigned (if an CNV exists at this locus);
#          *%maxP* - the confidence of mutation assignment.
#          *f* - Deprecated. The maximum likelyhood cellular prevalence
#          of this point mutation, before it has been assigned to SP.
#          This value is based on the copy number and allele frequency
#          of the mutation exclusively and is independent of other point
#          mutations. Column SP is less sensitive to noise and
#          considered the more accurate estimation of cellular mutation
#          prevalence.
#          *PM* - the total ploidy of all alleles in the subpopulation
#          harboring the point mutation (SP).
#          *PM_B* - the ploidy of the B-allele in the subpopulation
#          harboring the point mutation (SP).
#          *PM_cnv* - the total ploidy of all alleles in the
#          subpopulation harboring an CNV (SP_cnv).
#          *PM_B_cnv* - the ploidy of the B-allele, in the CNV harboring
#          subpopulation (SP_cnv).
#          If phylogeny reconstruction was successful, matrix includes
#          one additional column for each subpopulation from the
#          phylogeny, indicating whether or not the point mutation is
#          present in the corresponding subpopulation.

#densities: Matrix as obtained by 'computeCellFrequencyDistributions.'
#          Each row corresponds to a mutation and each column
#          corresponds to a cellular frequency. Each value
#          densities[i,j] represents the probability that mutation i is
#          present in a fraction f of cells, where f is given by:
#          colnames(densities[,j]).

#  ploidy: Matrix as obtained by 'assignQuantityToSP.' Each row
#          corresponds to a copy number segment, e.g. as obtained from a
#          circular binary segmentation algorithm. Includes one
#          additional column for each predicted subpopulation,
#          containing the ploidy of each segment in the corresponding
#          subpopulation.

#    tree: An object of class "phylo" (library ape) as obtained by
#          'buildPhylo.' Contains the inferred phylogenetic
#          relationships between subpopulations.

#结果文件
# 1. out.expands.spstats  包含了亚克隆的信息
#   Mean Weighted   score   precision       nMutations      Ancestor        ClosestDescendant
#   0.347571968513367       0.450226170010756       0.0412619947522278      36      0.636405931778961       NA
#   0.636405931778961       1.65155788590212        0.0412619947522278      20      0.8427159055401 0.347571968513367
#   0.8427159055401 1.34021083604764        0.0412619947522278      32      NA      0.636405931778961

# 2. displaying a visual representation of the identified subpopulations (1. subpopulation size   2. phylogeny)



##############多点取样####################
# Inferring phylogenetic relations between subpopulations from multiple geographical tumor samples
# 需要输入文件：
#   The CBS files for each sample: out.expands.sps.cbs
#   The SP files for each sample (previously calculated via runExPANdS-function): out.expands.sps