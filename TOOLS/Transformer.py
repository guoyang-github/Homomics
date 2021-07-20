# coding=utf-8
# guoyang

import argparse
import os.path
import sys
import HTSeq
import collections
import math
import Vcf
import re

usage='''
various types of format transfor
'''


def annoVcf2MafForAbsolute(args):
	'''
	适用于absolute的maf生成有两种方式：
		1. 用vcf2maf工具，结果直接可用，不过太重量级  
		2. 从annovar注释过（gene-based annotatition first）的VCF中用此函数转换
	'''
	Iterms = ['Hugo_Symbol', 'Chromosome', 'Start_position', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 't_ref_count', 't_alt_count']
	result = []
	vcfRecord = Vcf.VcfReader(args.annoVcf)
	for eachRecord in vcfRecord:
		AD = eachRecord.SAMPLE[args.tumorId]['AD'].split(',')   #Allelic depths for the ref and alt alleles in the order listed
		tempRecord = ['Gene.' + args.table, eachRecord.CHROM, eachRecord.POS, '', args.tumorId, AD[0], AD[1]]
		result.append('\t'.join(tempRecord) + '\n')
	outfile = os.path.join(args.outdir, os.path.basename(args.annoVcf).replace('.vcf', '') + '.maf')
	with open(outfile, 'w') as outFile:
		outFile.writelines(['\t'.join(Iterms) + '\n'] + result)

def annoVcf2MultiannoTxtForMafTools(args):
	'''
	适用于maftools的multianno.txt生成有两种方式：
		1. 用vcf2maf工具，结果直接生成maf可用，不过太重量级  
		2. 从annovar注释过（gene-based annotatition first）的VCF中用此函数转换出来multianno.txt,再用annovar.maf函数，得到maftools适用的maf，支持多样本合并
	'''
	Iterms = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene.' + args.table, 'GeneDetail.' + args.table, 'ExonicFunc.' + args.table, 'AAChange.' + args.table, 'Tumor_Sample_Barcode', 'Func.' + args.table]
	sampleVcf = {}
	with open(args.annoVcf) as sampleFile:
		for eachLine in sampleFile:
			if eachLine.strip() == '':
				continue
			temp = eachLine.strip().split()
			sampleVcf[temp[0]] = temp[1]
	result = []
	for eachSample in sampleVcf:
		vcfRecord = Vcf.VcfReader(sampleVcf[eachSample])
		for eachRecord in vcfRecord:
			tempRecord = [eachRecord.CHROM, eachRecord.start, eachRecord.end, eachRecord.REF, eachRecord.ALT, eachRecord.INFO['Gene.' + args.table], eachRecord.INFO['GeneDetail.' + args.table], eachRecord.INFO['ExonicFunc.' + args.table], eachRecord.INFO['AAChange.' + args.table], eachSample, eachRecord.INFO['Func.' + args.table]]
			result.append('\t'.join(tempRecord) + '\n')
	outfile = os.path.join(args.outdir, os.path.basename(args.annoVcf) + '.multianno.txt')
	with open(outfile, 'w') as outFile:
		outFile.writelines(['\t'.join(Iterm)+ '\n' ] + result)

def annoVcf2MultiannoTxtForOncodriveCLUST(args):
	'''
	annovar的注释结果中，突变位点会包含在同一个基因的不同转录本上，所以这个选择就有些问题
	本函数的策略是选择AAchange中一串描述中，p.X???Y 中，蛋白质序列???最长的那个
	'''
	Iterms = ['GeneSymbol', 'Transcript', 'Sample', 'ExonicFunc', 'AAchange', 'ProteinPosition']
	sampleVcf = {}
	with open(args.annoVcf) as sampleFile:
		for eachLine in sampleFile:
			if eachLine.strip() == '':
				continue
			temp = eachLine.strip().split()
			sampleVcf[temp[0]] = temp[1]

	def getAAchange(aachange):
		aachanges = aachange.strip().split(',')
		maxPosition = 0
		maxAAChange = ''
		AAChange = {}
		for eachAAChange in aachanges:
			if eachAAChange not in AAChange:
				AAChange[eachAAChange] = {}
			AAChange[eachAAChange]['Transcript'] = eachAAChange.split(':')[1]
			AAChange[eachAAChange]['proteinPosition'] = re.search(r'p\.\w(\d+)\w', eachAAChange).group(1)
			if int(AAChange[eachAAChange]['proteinPosition']) >= maxPosition:
				maxPosition = int(AAChange[eachAAChange]['proteinPosition'])
				maxAAChange = eachAAChange
		return AAChange[maxAAChange]['Transcript'], maxAAChange, AAChange[maxAAChange]['proteinPosition']    #返回蛋白质序列最大的AAChange

	resultSyn = []
	resultNonSyn = []
	for eachSample in sampleVcf:
		vcfRecord = Vcf.VcfReader(sampleVcf[eachSample])
		for eachRecord in vcfRecord:
			if eachRecord.INFO['Func.' + args.table] != 'exonic' or eachRecord.INFO['ExonicFunc.' + args.table] != 'unknown':
				continue
			genesymbol = eachRecord.INFO['Gene.' + args.table].strip().split(',')[0]
			transcript, aachange, proteinposition = getAAchange(eachRecord.INFO['AAChange.' + args.table])
			exonicfunc = eachRecord.INFO['ExonicFunc.' + args.table]
			tempRecord = [genesymbol, transcript, eachSample, exonicfunc, aachange, proteinposition]
			if exonicfunc.startswith('synonymous'):
				resultSyn.append('\t'.join(tempRecord) + '\n')
			else:
				resultNonSyn.append('\t'.join(tempRecord) + '\n')
	
	outsyn = os.path.join(args.outdir, os.path.basename(args.annoVcf) + '.syn.txt')
	outnonsyn = os.path.join(args.outdir, os.path.basename(args.annoVcf) + '.nonsyn.txt')
	with open(outsyn, 'w') as outSyn:
		outSyn.writelines(resultSyn)
	with open(outnonsyn, 'w') as outNonSyn:
		outNonSyn.writelines(resultNonSyn)			

def vepTab2MultiannoTxtForOncodriveCLUST(args):
	'''
	从VEP的注释结果可以更可靠的得到结果
	VEP采用--pick 只保留canonical的转录本结果，注释结果里包括CDS长度／蛋白质位置
	--pick_order canonical,appris,tsl,biotype,ccds,rank,length
	VEP --tab 的表格结果，很直观
	'''
	Iterms = ['GeneSymbol', 'Gene', 'Transcript', 'Sample', 'Consequence', 'ProteinPosition']
	ConsequenceIterms = ['upstream_gene_variant', 'intron_variant', 'synonymous_variant', 'missense_variant', 'splice_region_variant', '3_prime_UTR_variant', 'non_coding_transcript_variant', 'stop_gained', 'NMD_transcript_variant', '5_prime_UTR_variant', 'downstream_gene_variant', 'intergenic_variant', 'non_coding_transcript_exon_variant', 'stop_lost', 'stop_retained_variant', 'start_lost', 'regulatory_region_variant']
	exonicConsequenceIterms = ['synonymous_variant', 'missense_variant', 'stop_gained', 'stop_lost']
	synConsequenceIterms = ['synonymous_variant']
	nonsynConsequenceIterms = ['missense_variant', 'stop_gained', 'stop_lost']
	gtIterms = ['Symbol', 'Ensembl Transcript ID', 'CDS Length']
	sampleTab = {}
	with open(args.vepTab) as sampleFile:
		for eachLine in sampleFile:
			if eachLine.strip() == '':
				continue
			temp = eachLine.strip().split()
			sampleTab[temp[0]] = temp[1]
	resultSyn = []
	resultNonSyn = []
	resultGT = {}
	for eachSample in sampleTab:
		with open(sampleTab[eachSample]) as tabFile:
			for eachLine in tabFile:
				if eachLine.strip() == '' or eachLine.startswith('##'):
					continue
				if 	eachLine.startswith('#Uploaded_variation'):
					temp = eachLine.strip().split('\t')
					GeneIndex = temp.index('Gene')
					FeatureIndex = temp.index('Feature')
					ConsequenceIndex = temp.index('Consequence')
					CDS_positionIndex = temp.index('CDS_position')
					Protein_positionIndex = temp.index('Protein_position')
					SYMBOLIndex = temp.index('SYMBOL')
					CANONICALIndex = temp.index('CANONICAL')
					continue
				temp = eachLine.strip().split('\t')
				Consequence = temp[ConsequenceIndex].split(',')
				if len(set(Consequence) & set(exonicConsequenceIterms)) == 0:
					continue
				GeneSymbol =  temp[SYMBOLIndex]
				Gene = temp[GeneIndex]
				Transcript = temp[FeatureIndex]
				ProteinPosition = temp[Protein_positionIndex].split('/')[0]
				CDSLength = temp[CDS_positionIndex].split('/')[1]
				CANONICAL = temp[CANONICALIndex]
				tempResult = [GeneSymbol, Gene, Transcript, eachSample, Consequence, ProteinPosition]
				if len(set(Consequence) & set(synConsequenceIterms)):
					resultSyn.append('\t'.join(tempResult) + '\n')
				else:
					resultNonSyn.append('\t'.join(tempResult) + '\n')
				if GeneSymbol not in resultGT:   #收集canonical，不是的位点也要记录，没有办法，但毕竟是少数
					resultGT[GeneSymbol] = {}
					resultGT[GeneSymbol]['Transcript'] = Transcript
					resultGT[GeneSymbol]['CDSLength'] = CDSLength
					resultGT[GeneSymbol]['CANONICAL'] = CANONICAL
				else:   #发现canonical的就替换原来不是canonical的
					if resultGT[GeneSymbol]['CANONICAL'] == 'NO' and CANONICAL == 'YES':
						resultGT[GeneSymbol]['Transcript'] = Transcript
						resultGT[GeneSymbol]['CDSLength'] = CDSLength
						resultGT[GeneSymbol]['CANONICAL'] = CANONICAL
	outsyn = os.path.join(args.outdir, os.path.basename(args.vepTab) + '.syn.txt')
	outnonsyn = os.path.join(args.outdir, os.path.basename(args.vepTab) + '.nonsyn.txt')   
	outgt = os.path.join(args.outdir, os.path.basename(args.vepTab) + '.gene_transcript.txt')
	with open(outsyn, 'w') as outSyn:
		outSyn.writelines(resultSyn)
	with open(outnonsyn, 'w') as outNonSyn:
		outNonSyn.writelines(resultNonSyn)
	with open(outgt, 'w') as outGT:
		result = []
		for eachGT in resultGT:
			result.append(eachGT + '\t' + resultGT[eachGT]['Transcript'] + '\t' + resultGT[eachGT]['CDSLength'] + '\n')	
		outGT.writelines(result)

def ControlFreeC4GenePattern(cnv, ratio, sample = None):
	'''
	利用controlfreec的CNV和ratio结果构建absolute/gistic2需要的segmentation输入文件
	(1)  Sample(sample name) (2)  Chromosome(chromosome number) (3)  Start Position (segment start position, in bases) (4)  End Position(segment end position, in bases) (5)  Num markers(number of markers in segment) (6)  Seg.CN(log2() -1 of copy number)
	median是相对于mean更加稳健的统计量，同时controlfreec提供的也是median，起初计划选用median，但发现同一个CNV可能会对应几个median值，也许是进行了多次的区间合并，所以按照absolute的要求，最后采用了mean值
	'''
	result = []
	with open(args.cnv) as regionSet:
		for eachRegion in regionSet:
			if eachRegion.strip() == '' or eachRegion.startswith('Chromosome'):
				continue
			probeCount = 0
			rawRatios = []
			#ratios = []   #收集median Ratio
			tempRegion =  eachRegion.strip().split()
			region = HTSeq.GenomicInterval(tempRegion[0], int(tempRegion[1])-1, int(tempRegion[2]))
			with open(args.ratio) as exomeSet:
				for eachExome in exomeSet:
					if eachExome.strip() == '' or eachExome.startswith('Chromosome'):
						continue
					tempExome = eachExome.strip().split()
					point = HTSeq.GenomicPosition(tempExome[0], int(tempExome[1])-1)
					if region.contains(point):
						probeCount += 1
						#ratios.append(float(tempExome[3]))
						rawRatios.append(float(tempExome[2]))
			#medianRatio, commonCount = collections.Counter(ratios).most_common(1)[0]
			meanRatio = sum(rawRatios)/len(rawRatios)
			if meanRatio == 0:   #if medianRatio == 0:   #发现数据文件_CNV.txt medianRatio列中有0值，放弃这些CNV
				continue
			#if float(commonCount)/probeCount < 0.9:
			#	sys.exit('different medianRatios in the same CNV region! \n' + eachRegion)
			tempResult = tempRegion[0] + '\t' + tempRegion[1] + '\t' + tempRegion[2] + '\t' + str(probeCount) + '\t' + str(math.log(float(meanRatio),2)) + '\n'
			if sample:
				result.append(sample + '\t' + tempResult)
			else:
				result.append(tempResult)
	return result

def ControlFreeC4Absolute(args):
	result = ControlFreeC4GenePattern(args.cnv, args.ratio)		
	outfile =  os.path.join(args.outdir, os.path.basename(args.cnv) + '.seg.dat')
	with open(outfile, 'w') as outFile:
		outFile.writelines(result)


def ControlFreeC4GISTIC2(args):
	segResult = []
	markerPervious =[]
	markerNext = []
	def makeMarkerFile(ratio):   #从_ratio.txt文件中提取 MarkerName Chromosome MarkerPosition
		markerResult = []
		with open(ratio) as ratioFile:
			for eachRatio in ratioFile:
				if eachRatio.strip() == '':
					continue
				tempRatio = eachRatio.strip().split()
				markerResult.append(tempRatio[0] + ':' + tempRatio[1] + '\t' + tempRatio[0] + '\t' + tempRatio[1] + '\n')
		return markerResult
	with open(args.sampleCnvRatio) as sampleCnvRatioFile:		
		for eachSample in sampleCnvRatioFile:
			if eachSample.strip() == '':
				continue
			tempSample = eachSample.strip().split()
			segResult += ControlFreeC4GenePattern(tempSample[1], tempSample[2], tempSample[0])
			if markerPervious:   #查查每个marker file是不是俩俩一致，只保留最后一个作为输出，节省内存，每个样本的都应该是一样的
				markerNext = makeMarkerFile(tempSample[2])
				if markerPervious == markerNext:
					markerPervious = markerNext
				else:
					sys.exit('ERROR: marker file is different between _ratio.txts!')
			else:
				markerPervious = makeMarkerFile(tempSample[2])
	segfile = os.path.join(args.outdir, os.path.basename(args.sampleCnvRatio) + '.seg.dat')
	markerfile = os.path.join(args.outdir, os.path.basename(args.sampleCnvRatio) + '.marker.dat')
	with open(segfile, 'w') as segFile:
		segFile.writelines(segResult)
	with open(markerfile, 'w') as markerFile:
		markerFile.writelines(markerNext)


def mutectControlFreeC4PyClone(args):
	Items = ['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn']
	CNVs={}   #获得CNV文件中的信息，包括total及allelic
	with open(args.cnv) as cnvFile:
		for eachCNV in cnvFile:
			if eachCNV.strip() == '':
				continue
			tempCNV = eachCNV.strip().split()
			cnvID = '_'.join(tempCNV[:3])
			if cnvID not in CNVs:
				CNVs[cnvID] = {}
			CNVs[cnvID]['region'] = HTSeq.GenomicInterval(tempCNV[0], int(tempCNV[1])-1, int(tempCNV[2]))
			CNVs[cnvID]['cn'] = tempCNV[3]
			#GT = collections.Counter(tempCNV[5]).most_common(2)   #2倍区域（不在CNV文件里的区域），基因型无法判断，其他区域也无法判断准确，所以还是不用allelic的CNV了，在pyClone.py脚本里具体提到了这个问题
			#CNVs[cnvID]['major_cn'] = GT[0][1]
			#CNVs[cnvID]['minor_cn'] = GT[1][1]
			#assert int(CNVs[cnvID]['cn']) == int(CNVs[cnvID]['major_cn']) + int(CNVs[cnvID]['minor_cn'])
	result = []
	vcfRecord = Vcf.VcfReader(args.vcf)
	for eachRecord in vcfRecord:
		major_cn = 2
		tempPosition = HTSeq.GenomicPosition(eachRecord.CHROM, int(eachRecord/POS)-1)
		for eachCNV in CNVs:
			if CNVs[eachCNV]['region'].contains(tempPosition):
				major_cn = CNVs[eachCNV]['cn']
				break
		AD = eachRecord.SAMPLE[args.tumorId]['AD'].split(',')   #Allelic depths for the ref and alt alleles in the order listed
		result.append(eachRecord.CHROM + '_' + eachRecord.POS + '\t' + AD[0] + '\t' + AD[1] + '\t' + '2' + '\t' + '0' + '\t' + str(major_cn) + '\n')
	outfile = os.path.join(args.outdir, os.path.basename(args.vcf) + '.pyclone.dat')
	with open(outfile, 'w') as outFile:
		outFile.writelines(['\t'.join(Iterm)+ '\n' ] + result)
	

def varscan4EXPANDS(args):
	#*REF* - ASCII code of the reference nucleotide (in hg18/hg19)
	#*ALT* - ASCII code of the B-allele nucleotide    还有一个信息从这里得到，所谓B-allele就是指和ref不一样的变异ALT
	#*AF_Tumor* - allele frequeny of B-allele
	#*PN_B* - ploidy of B-allele in normal cells. 
		#A value of 0 indicates that the mutation has only been detected in the tumor sample (i.e. somatic mutations). 
		#A value of 1 indicates that the mutation is also present in the normal (control) sample, albeit at reduced allele frequency (i.e. mutation is consequence of LOH). 
		#Other mutations should not be included.
	Items = ['chr', 'startpos', 'endpos', 'REF', 'ALT', 'AF_Tumor', 'PN_B']
	result = []
	def extractVcf(vcfFile, type):
		simVcf = []
		vcfRecord = Vcf.VcfReader(vcfFile)
		for eachRecord in vcfRecord:
			if float(eachRecord.SAMPLE['TUMOR']['FREQ'].strip('%')) > float(eachRecord.SAMPLE['NORMAL']['FREQ'].strip('%')):
				tempRecord = [eachRecord.CHROM, eachRecord.POS, eachRecord.POS, ord(eachRecord.REF), ord(eachRecord.ALT), eachRecord.SAMPLE['TUMOR']['FREQ'], type]
				simVcf.append('\t'.join(tempRecord) + '\n')
		return simVcf
	result += extractVcf(args.somatic, 0)
	if args.loh:
		result += extractVcf(args.loh, 1)
	outfile = os.path.join(args.outdir, 'varscan2.expands.snv')
	with open(outfile, 'w') as outFile:
		outFile.writelines(['\t'.join(Iterm)+ '\n' ] + result)


def controlfreec4EXPANDS(args):
	# controlfreec 和 exomeCNV的 CN值都是估计成整数的值
	# *CN_Estimate* - average copy-number of the segment among all cells. (absolute copy number)
	# 从expands的测试数据来看，cbs里面的CN_Estimate最好能是用copy number ratio乘以ploidy得到
	# controlfreec 可以参考函数ControlFreeC4GenePattern从_CNV _ratio两个文件中得到
	# 对于WES的情况，exomeCNV的结果文件里包括CN 和 copy number ratio 的值，可以直接取得
	Item = ['chr', 'startpos', 'endpos', 'CN_Estimate']
	result = []
	with open(args.cnv, 'w') as cnvFile:
		for eachCNV in cnvFile:
			tempItem = eachCNV.strip().split()
			result.append('\t'.join(tempItem[:4]) + '\n')
	outfile = os.path.join(args.outdir, 'controlfreec.expands.cbs')
	with open(outfile, 'w') as outFile:
		outFile.writelines(['\t'.join(Iterm)+ '\n' ] + result)


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	subparsers = parser.add_subparsers(help='various types of format transfor')
	#controlfreec4absolute
	parser_controlfreec4absolute = subparsers.add_parser('controlfreec4absolute', help='from controlfreec to absolute Segmented copy ratios data')
	parser_controlfreec4absolute.add_argument('--cnv', help = '_CNV', required = True)
	parser_controlfreec4absolute.add_argument('--ratio', help =' _ratio.txt', required = True)
	parser_controlfreec4absolute.add_argument('--outdir', help = 'output dir', default = '.')
	parser_controlfreec4absolute.set_defaults(func = ControlFreeC4Absolute)

	#ControlFreeC4GISTIC2
	parser_controlfreec4gistic2 = subparsers.add_parser('controlfreec4gistic2', help='from controlfreec to gistic2 Segmented copy ratios data')
	parser_controlfreec4gistic2.add_argument('--sampleCnvRatio', help = 'sampleName\t_CNV.txt\t_ratio.txt\n', required = True)
	parser_controlfreec4gistic2.add_argument('--outdir', help = 'output dir', default = '.')
	parser_controlfreec4gistic2.set_defaults(func = ControlFreeC4GISTIC2)

	#mutectcontrolfreec4pyclone
	parser_mutectcontrolfreec4pyclone = subparsers.add_parser('mutectcontrolfreec4pyclone', help='from mutect and controlfreec to pyclone')
	parser_mutectcontrolfreec4pyclone.add_argument('--cnv', help = 'controlfreec _CNV', required = True)
	parser_mutectcontrolfreec4pyclone.add_argument('--vcf', help ='mutect vcf', required = True)
	parser_mutectcontrolfreec4pyclone.add_argument('--tumorId', help ='tumor ID', required = True)
	#parser_mutectcontrolfreec4pyclone.add_argument('--type', help = 'toal CNV or allelic CNV', choices = ['total', 'allelic'], default = 'total')
	parser_mutectcontrolfreec4pyclone.add_argument('--outdir', help = 'output dir', default = '.')
	parser_mutectcontrolfreec4pyclone.set_defaults(func = mutectControlFreeC4PyClone)	

	#annoVcf2MafForAbsolute
	parser_annoVcf2MafForAbsolute = subparsers.add_parser('annoVcf2MafForAbsolute', help='from annotated VCF to Maf for ABSOLUTE')
	parser_annoVcf2MafForAbsolute.add_argument('--annoVcf', help = 'annotated VCF', required = True)
	parser_annoVcf2MafForAbsolute.add_argument('--tumorId', help = 'tumor ID', required = True)
	parser_annoVcf2MafForAbsolute.add_argument('--outdir', help = 'output dir', default = '.')
	parser_annoVcf2MafForAbsolute.add_argument('--table', help = "reference table used for gene-based annotations. Can be 'ensGene' or 'refGene'. Default 'refGene'", choices = ['refGene', 'ensGene'], default = 'refGene')
	parser_annoVcf2MafForAbsolute.set_defaults(func = annoVcf2MafForAbsolute)	

	#annoVcf2MultiannoTxtForMafTools
	parser_annoVcf2MultiannoTxtForMafTools = subparsers.add_parser('annoVcf2MultiannoTxtForMafTools', help='from annotated VCF to Maf for maftools')
	parser_annoVcf2MultiannoTxtForMafTools.add_argument('--annoVcf', help = "tumorID\tannoVcf\ntumorID\tannoVcf\n", required = True)
	parser_annoVcf2MultiannoTxtForMafTools.add_argument('--outdir', help = "output dir", default = '.')
	parser_annoVcf2MultiannoTxtForMafTools.add_argument('--table', help = "reference table used for gene-based annotations. Can be 'ensGene' or 'refGene'. Default 'refGene'", choices = ['refGene', 'ensGene'], default = 'refGene')
	parser_annoVcf2MultiannoTxtForMafTools.set_defaults(func = annoVcf2MultiannoTxtForMafTools)	

	#annoVcf2MultiannoTxtForOncodriveCLUST
	parser_annoVcf2MultiannoTxtForOncodriveCLUST = subparsers.add_parser('annoVcf2MultiannoTxtForOncodriveCLUST', help='from mutect annotated VCF to txt for OncodriveCLUST')
	parser_annoVcf2MultiannoTxtForOncodriveCLUST.add_argument('--annoVcf', help = "tumorID\tannoVcf\ntumorID\tannoVcf\n", required = True)
	parser_annoVcf2MultiannoTxtForOncodriveCLUST.add_argument('--outdir', help = "output dir", default = '.')
	parser_annoVcf2MultiannoTxtForMafTools.add_argument('--table', help = "reference table used for gene-based annotations. Can be 'ensGene' or 'refGene'. Default 'refGene'", choices = ['refGene', 'ensGene'], default = 'refGene')
	parser_annoVcf2MultiannoTxtForOncodriveCLUST.set_defaults(func = annoVcf2MultiannoTxtForOncodriveCLUST)

	#vepTab2MultiannoTxtForOncodriveCLUST
	parser_vepTab2MultiannoTxtForOncodriveCLUST = subparsers.add_parser('vepTab2MultiannoTxtForOncodriveCLUST', help='from mutect annotated VCF to txt for OncodriveCLUST')
	parser_vepTab2MultiannoTxtForOncodriveCLUST.add_argument('--vepTab', help = "tumorID\tvepTab\ntumorID\tvepTab\n", required = True)
	parser_vepTab2MultiannoTxtForOncodriveCLUST.add_argument('--outdir', help = "output dir", default = '.')
	parser_vepTab2MultiannoTxtForOncodriveCLUST.set_defaults(func = vepTab2MultiannoTxtForOncodriveCLUST)

	#varscan4EXPANDS
	parser_varscan4EXPANDS = subparsers.add_parser('varscan4EXPANDS', help='from varscan2 somatic and LOH point mutation file to snv for EXPANDS')
	parser_varscan4EXPANDS.add_argument('--somatic', help = 'snpSomaticFilterfpFilter', required = True)
	parser_varscan4EXPANDS.add_argument('--loh', help = 'snpLOHfpFilter', default = None)
	parser_varscan4EXPANDS.add_argument('--outdir', help = 'output dir', default = '.')

	#controlfreec4EXPANDS
	parser_controlfreec4EXPANDS = subparsers.add_parser('controlfreec4EXPANDS', help='from controlfreec CN to cbs for EXPANDS')
	parser_controlfreec4EXPANDS.add_argument('--cnv', help = 'controlfreec _CNV', required = True)
	parser_controlfreec4EXPANDS.add_argument('--outdir', help = 'output dir', default = '.')

	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	args.func(args)
	