# coding=utf-8
# guoyang
import sys 
import os.path
import re

class VcfRecord:
	def __init__(self,record,samples):
		self.record = record
		self.samples = samples
		self.Record = self.record.strip().split()
		self.CHROM = self.Record[0]
		self.POS = int(self.Record[1])
		self.ID = self.Record[2]
		self.REF = self.Record[3]
		self.ALT = self.Record[4].split(',')
		self.QUAL = self.Record[5]
		self.FILTER = self.Record[6]
		self.INFO = self._getInfo()
		self.FORMAT = self.Record[8].split(':')
		self.SAMPLE = self._getSample()
		self.start = self.POS   #Maf format 1-based
		self.end = self.POS + len(self.REF) - 1  #Maf format 1-based end inclusive
		self.type = self._getType()

	#Get each Item in INFO field   Value: dict
	def _getInfo(self):
		tempInfo = {}
		for eachTerm in self.Record[7].split(';'):
			if '=' not in eachTerm:
				continue
			tempeachTerm = eachTerm.split('=')
			tempInfo[tempeachTerm[0]] = tempeachTerm[-1]		
		return tempInfo
	
	#Get each Item in FORMAT field for each sample   Value: dict
	def _getSample(self):
		tempSample = {}
		if len(self.samples) != len(self.Record[9:]):
			sys.exit("ERROR: NOT EQUAL -- number of samples vs FORMAT!")
		for i, eachCall in enumerate(self.Record[9:]):
			if self.samples[i] not in tempSample:
				tempSample[self.samples[i]] = {}
			else:
				sys.exit('ERROR: SAME sample name for different FORMAT!')
			callTerms = eachCall.split(':')
			if len(callTerms) != len(self.FORMAT):
				sys.exit('ERROR: FORMAT terms deficiency')
			for j,eachCallTerm in enumerate(callTerms):
				if self.FORMAT[j] not in tempSample[self.samples[i]]:
					tempSample[self.samples[i]][self.FORMAT[j]] = eachCallTerm
				else:
					sys.exit('ERROR: duplicated FORMAT term!')	
		return tempSample

	def _getType(self):
		if len(self.REF) == 1:
			for eachAlt in self.ALT:
				if len(eachAlt) > 1:
					return 'INS'
			return 'SNP'
		else:
			for eachAlt in self.ALT:
				if len(eachAlt) > len(self.REF):
					return 'INS'
				if len(eachAlt) < len(self.REF):
					return 'DEL'
	
	def genotype(self, sample):
		if  sample not in self.samples:
			sys.exit('ERROR: illegal sample name!')
		if 'GT' not in self.SAMPLE[sample]:
			return None
		GT = re.split(r'[/|]', self.SAMPLE[sample]['GT'])
		numGT = map(int,GT)
		if len(self.ALT) != max(numGT):
			sys.exit('ERROR: illegal genotype number!')
		genotypeSet = [self.REF]+self.ALT
		genoType = []
		for eachGT in numGT:
			genoType.append(genotypeSet[eachGT])
		return '/'.join(genoType)



#####################################################################	
class VcfReader:
	def __init__(self,vcf):
		self.vcf = vcf
		self.fileHandle = self._safeOpen() #Get the file object safely
		self.fileformat = ''
		self.metaInfo = {}
		self.reference = ''
		self.samples = []
		self._getMeta()
	
	#destructor to close the file object
	def __del__(self):
		self.fileHandle.close()
	
	#Due to the existence of comma, use this function to sperate the information in <>
	def _specialSplit(self, data, sign):
		last = ''   
		metas = []
		for eachItem in data.split(sign):
			if '=' in eachItem:
				metas.append(last)
				last = eachItem
			else:
				last += ' ,' + eachItem
		metas.append(last)
		metas.pop(0)
		return metas		


	#Get meta information, and turn the pointer to main body
	def _getMeta(self):
		for eachLine in self.fileHandle:
			if not eachLine.strip():
				continue
			if eachLine.startswith('##fileformat'):
				self.fileformat = eachLine.replace('##fileformat=', '')
			pattern = re.search(r'##(.+?)=<(.+)>', eachLine)
			if pattern:
				if pattern.group(1) not in self.metaInfo:
					self.metaInfo[pattern.group(1)] = {}

				metas = self._specialSplit(pattern.group(2), ',')
				
				for eachElement in metas:
					elements = eachElement.split('=', 1)
					if 'ID' in elements and elements[1] not in self.metaInfo[pattern.group(1)]:
						ID = elements[1]
						self.metaInfo[pattern.group(1)][ID] = {}
					self.metaInfo[pattern.group(1)][ID][elements[0]] = elements[1]

			if eachLine.startswith('##reference'):
				self.reference = eachLine.replace('##reference=file://', '')
			
			if eachLine.startswith('#CHROM'):
				temp = eachLine.strip().split()
				self.samples = temp[temp.index('FORMAT') + 1:] #Get sample list
				break
	
	def __iter__(self):
		return self
	
	def next(self):
		return VcfRecord(self.fileHandle.next(),self.samples)
	
	#Private method to open the VCF file, either zipped or unzipped
	def _safeOpen(self):
		try:
			if self.vcf.endswith('.gz'):
				import gzip
				return gzip.open(self.vcf)
			else:
				return open(self.vcf)
		except:
			sys.exit('ERROR: file not found!')



if __name__ == '__main__':
	vcf = 'test.vcf'
	vcfRecord = Vcf.VcfReader(vcf)
	print 'reference: ' + vcfRecord.reference
	print 'Depth: ' + vcfRecord.metaInfo['FORMAT']['DP']['Description']
	print 'X length: ' + vcfRecord.metaInfo['contig']['X']['length']
	print 'samples: ' + ' '.join(vcfRecord.samples)
	for eachRecord in vcfRecord:
		print 'CHROM: ' + eachRecord.CHROM
		print 'POS: ' + str(eachRecord.POS)
		print 'ALT: ' + ','.join(eachRecord.ALT)
		print 'GT: ' + eachRecord.SAMPLE['H625_S1E']['GT']
		break