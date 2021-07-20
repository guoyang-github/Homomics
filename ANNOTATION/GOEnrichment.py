from scipy import stats
import argparse
#guoyang 2014.12.18

parser = argparse.ArgumentParser(description = "GO enrichment guoyang@novogene.com")
parser.add_argument('--genelist', help = "gene list for GO enrichment")
parser.add_argument('--go', help = "info from ensembl biomart with format: gene\tGO term\tGO description")
parser.add_argument('--cutoff', help = "significant level default: 0.05", default = 0.05, type = float)

argv = parser.parse_args()

geneList = set()
for eachLine in open(argv.genelist):
    if eachLine.strip() == '':
        continue
    geneList.add(eachLine.strip())

go_file = open(argv.go)
go_file.readline()    
goMapgene = {}
geneMapgo = {}
goDes = {}
backgroud = set()
for eachLine in go_file:
    if eachLine.strip() == '':
        continue
    temp=eachLine.strip().split('\t')
    if len(temp) == 1:
        continue
    backgroud.add(temp[0])
    if temp[1] not in goDes:
        goDes[temp[1]]=temp[2]
    if temp[0] not in geneMapgo:
        geneMapgo[temp[0]] = set()
    geneMapgo[temp[0]].add(temp[1])
    if temp[1] not in goMapgene:
        goMapgene[temp[1]] = set()
    goMapgene[temp[1]].add(temp[0])
go_file.close()

Map=[]
for eachGene in geneMapgo:
    Map.append(eachGene + '\t' + '\t'.join(geneMapgo[eachGene]) + '\n')
open('gene.go', 'w').writelines(Map)

Map=[]
for eachGO in goMapgene:
    Map.append(eachGO + '\t' + '\t'.join(goMapgene[eachGO]) + '\n')
open('go.gene', 'w').writelines(Map)


    
print '%s of geneList in backgroud with GO annotation' % (float(len(geneList & backgroud)) / len(geneList))
if len(geneList & backgroud) < float(len(geneList))/2:
    print 'WARNING: more than half of genes in geneList are not in backgroud !'


result_all = {}
for eachGO in goMapgene:
    M = len(backgroud)
    N = len(geneList - backgroud)
    n = len(goMapgene[eachGO])
    hg = stats.hypergeom(M, n, N)
    k = len(geneList & goMapgene[eachGO])
    p = 1 - hg.cdf(k)
    if p not in result_all:
        result_all[p] = []
    result_all[p].append(eachGO + '\t' + goDes[eachGO] + '\t' + \
                        ','.join(geneList & goMapgene[eachGO]) + \
                        '\t' + str(p) + '\n')

result = []
result_sig = []
for eachP in sorted(result_all.keys()):
    result += result_all[eachP]
    if eachP <= argv.cutoff:
        result_sig += result_all[eachP]
    
open('GO_enrichment_list.txt', 'w').writelines(result)
open('Go_enrichment_sig.txt', 'w').writelines(result_sig)