import HTSeq
import os.path
import argparse
import itertools
import gzip


usage = '''
Demultiplexing for single index
'''

def paramsParse():
    parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--indexFile', help = "index file   samplename\tindexsequence\n", required = True)
    parser.add_argument('--read1', help = "Read1 file", required = True)
    parser.add_argument('--read2', help = "Read2 file", required = True)
    parser.add_argument('--outDir', help = "output dir", default ='.')
    parser.add_argument('--mismatch', help = "mismatch for index comparation", default = 1)
    return parser.parse_args()


def hamming_distance(s1, s2):
    "Return the Hamming distance between equal-length sequences"
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(el1 != el2 for el1, el2 in zip(s1, s2))


def index(indexFile):
    sampleIndex = {}
    for eachLine in open(indexFile):
        if eachLine.strip() == '':
            continue
        sample, index = eachLine.strip().split()
        if sample not in sampleIndex:
            sampleIndex[sample] = {}
        sampleIndex[sample]['index'] = index
    return sampleIndex


def demultiplexing(sampleIndex, read1, read2, outdir, mismatch):
    Read1 = HTSeq.FastqReader(read1)
    Read2 = HTSeq.FastqReader(read2)
    for eachSample in sampleIndex:
        sampleIndex[eachSample]['Read1'] = gzip.open(os.path.join(outdir, eachSample + '_' + sampleIndex[eachSample]['index'] + '_1.fq.gz'), mode = 'wb', compresslevel = 1)
        sampleIndex[eachSample]['Read2'] = gzip.open(os.path.join(outdir, eachSample + '_' + sampleIndex[eachSample]['index'] + '_2.fq.gz'), mode = 'wb', compresslevel = 1)
    undetermined1 = gzip.open(os.path.join(outdir, 'undetermined_1.fq.gz'), mode = 'wb', compresslevel = 1)
    undetermined2 = gzip.open(os.path.join(outdir, 'undetermined_2.fq.gz'), mode = 'wb', compresslevel = 1)
    for R1, R2 in itertools.izip(Read1, Read2):
        indexSeq = R1.name.split(':')[-1].strip()
        bestMatch = None
        bestMisCount = 1000
        for eachSample in sampleIndex:
            HD = hamming_distance(indexSeq, sampleIndex[eachSample]['index'])
            if  HD == 0:
                bestMatch = eachSample
                break
            if HD <= mismatch and HD <= bestMisCount:
                bestMatch = eachSample
                bestMisCount = HD
        if bestMatch:
            R1.write_to_fastq_file(sampleIndex[bestMatch]['Read1'])
            R2.write_to_fastq_file(sampleIndex[bestMatch]['Read2'])
        else:
            R1.write_to_fastq_file(undetermined1)
            R2.write_to_fastq_file(undetermined2)
    for eachSample in sampleIndex:
        sampleIndex[eachSample]['Read1'].close()
        sampleIndex[eachSample]['Read2'].close()
    undetermined1.close()
    undetermined2.close()


if __name__ == '__main__':
    args=paramsParse()
    demultiplexing(index(args.indexFile), args.read1, args.read2, args.outDir, args.mismatch)
