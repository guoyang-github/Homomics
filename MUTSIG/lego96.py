# coding=utf-8
# guoyang
import HTSeq
import argparse
import Vcf
import collections
import os.path

usage = '''
从多个样本的vcf中提取tri-nucleotide mutation contexts
'''


triNuclMutCont = ["TGTT", "TGGT", "TGCT", "TGAT", "TGTG", "TGGG", "TGCG", "TGAG", "TGTC", "TGGC", "TGCC", "TGAC", "TGTA", "TGGA", "TGCA", "TGAA", "TCTT", "TCGT", "TCCT", "TCAT", "TCTG", "TCGG", "TCCG", "TCAG", "TCTC", "TCGC", "TCCC", "TCAC", "TCTA", "TCGA", "TCCA", "TCAA", "TATT", "TAGT", "TACT", "TAAT", "TATG", "TAGG", "TACG", "TAAG", "TATC", "TAGC", "TACC", "TAAC", "TATA", "TAGA", "TACA", "TAAA", "CAAA", "CAAC", "CAAG", "CAAT", "CACA", "CACC", "CACG", "CACT", "CAGA", "CAGC", "CAGG", "CAGT", "CATA", "CATC", "CATG", "CATT", "CGAA", "CGAC", "CGAG", "CGAT", "CGCA", "CGCC", "CGCG", "CGCT", "CGGA", "CGGC", "CGGG", "CGGT", "CGTA", "CGTC", "CGTG", "CGTT", "CTAA", "CTAC", "CTAG", "CTAT", "CTCA", "CTCC", "CTCG", "CTCT", "CTGA", "CTGC", "CTGG", "CTGT", "CTTA", "CTTC", "CTTG", "CTTT"]


class Point:
    def __init__(self, pos, ref, alt, upstream = None, downstream = None):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.upstream = upstream
        self.downstream = downstream
    
    def getContext(self):
        if self.upstream and self.downstream:
            return (self.ref + self.alt + self.upstream + self.downstream).upper()
        else:
            return None


def extractVcf(vcf):
    vcfPoint = {}
    vcfRecord = Vcf.VcfReader(vcf)
    for eachReacord in vcfRecord:
        if eachReacord.CHROM not in vcfPoint:
            vcfPoint[eachReacord.CHROM] = []
        vcfPoint[eachReacord.CHROM].append(Point(eachReacord.POS, eachReacord.REF, eachReacord.ALT[0]))
    return vcfPoint


def extractPoint(nameVcf, reference):
    pointSet = collections.OrderedDict()
    with open(nameVcf) as mapfile:
        for eachLine in mapfile:
            if not eachLine:
                continue
            temp = eachLine.strip().split()
            sample = temp[0]
            vcf = temp[1]
            pointSet[sample] = extractVcf(vcf)
    for eachSeq in HTSeq.FastaReader(reference):
        for eachSample in pointSet:
            if eachSeq.name not in pointSet[eachSample]:
                continue
            for eachPoint in pointSet[eachSample][eachSeq.name]:
                assert eachPoint.ref == eachSeq.seq[eachPoint.pos-1]
                eachPoint.upstream = eachSeq.seq[eachPoint.pos-1-1]
                eachPoint.downstream = eachSeq.seq[eachPoint.pos-1+1]
    return pointSet


def statMatrix(pointSet):
    sMatrix = collections.OrderedDict()
    for eachSample in pointSet:
        sMatrix[eachSample] = {}
        for eachCont in triNuclMutCont:
            sMatrix[eachSample][eachCont] = 0
        for eachChrom in pointSet[eachSample]:
            for eachPoint in pointSet[eachSample][eachChrom]:
                if eachPoint.getContext() in triNuclMutCont:
                    sMatrix[eachSample][eachPoint.getContext()] += 1
    return sMatrix
            

def writeOutput(sMatrix, outdir):
    result = []
    sampleList = list(sMatrix)
    result.append('\t' + '\t'.join(sampleList) + '\n')
    for eachCont in triNuclMutCont:
        tempRow = []
        tempRow.append(eachCont)
        for eachSample in sampleList:
            tempRow.append(str(sMatrix[eachSample][eachCont]))
        result.append('\t'.join(tempRow) + '\n')
    with open(os.path.join(os.path.abspath(outdir), 'lego96.dat'), 'w') as output:
        output.writelines(result)


def paramsParse():
    parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--nameVcf', help = "Name\tSNPCVF\nName\tSNPCVF\n...",required = True)
    parser.add_argument('--reference', help = "genome reference file", required = True)
    parser.add_argument('--outdir', help = "output dir", default = '.')
    return parser.parse_args()



if __name__ == '__main__':
    args = paramsParse()
    pointSet = extractPoint(args.nameVcf, args.reference)
    sMatrix = statMatrix(pointSet)
    writeOutput(sMatrix, args.outdir)