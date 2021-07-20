# coding=utf-8
# guoyang
import subprocess
import argparse
import os.path
import sys
import os

usage='''
Somatic mutations are present in all cells of the human body and occur throughout life. 
They are the consequence of multiple mutational processes, including the intrinsic slight 
infidelity of the DNA replication machinery, exogenous or endogenous mutagen exposures, 
enzymatic modification of DNA and defective DNA repair. Different mutational processes 
generate unique combinations of mutation types, termed “Mutational Signatures”. 


Mutation Signature Profiling     http://archive.broadinstitute.org/cancer/cga/msp

What is a mutational process?
A catalogue of somatic mutations in cancer genomes is the outcome of cumulative actions of several mutagenic processes 
that operate over the patient’s lifetime, arising from an exposure to exogeneous DNA damaging agents or carcinogens 
(tobacco smoking / UV radiation / viral infections) or genomic defects in DNA repair or replicative processes 
(microsatellite instability / homologous recombination / DNA-polymerase proofreading defects) or endogeneous mutagens 
(reactive oxygen species or AID/APOBEC cytidine deaminases). [1,2] Even with a significant variation of mutation rates and 
a heterogeneity of mutation landscapes across tumor types somatic mutations often leave a characteristic molecular imprint, 
termed “mutation signatures”, i.e., a molecular profiling of mutational patterns in mutation types or local sequence contexts [3,4].

Bayesian non-negative matrix factorization (BayesNMF)
Characterizing underlying mutational processes with correct inferences for signature activities across samples provide
a key understanding on cancer initiation and progression. However, the number of mutation processes (K*) is highly variable 
across patients even in a single tumor type and its accurate estimation is a non-trivial task due to a different duration and 
intensity of exposure to a specific mutational process. Non-negative matrix factorization (NMF) has been widely used in deciphering 
mutations signatures in cancer somatic mutations stratified by 96 base substitutions in tri-nucleotide sequence contexts. In contrast 
to conventional NMF requiring K* as an input parameter BayesNMF [5,6,7] exploits the "automatic relevance determination" technique to 
infer the optimal K* from data itself at a balance between the “data fidelity” (likelihood) and the “model complexity” (regularization). 
As an example, Figure1A describes six mutational processes extracted from 14 WES mutation data of immortalized cell lines [8] and Figure1B 
represents the activity of discovered mutational processes across samples.
'''

class innerSoftware:
	def __init__(self):
		self.R_2_15_3 = 'Rscript'
		self.python = 'python'
		
class innerDatabase:
	def __init__(self):
		self.sigSet = ''
		self.reference = 'human_g1k_v37_decoy.fasta'	

class BAYESNMF:
	def __init__(self, nameVcf, outDir, tumorType, software = innerSoftware(), database = innerDatabase()):
		self.nameVcf = nameVcf
		self.tumorType = tumorType
		self.outDir = os.path.abspath(outDir)
		self.software = software
		self.database = database
		
	def LeGo96(self):
		self.software.lego96 = os.path.join(sys.path[0], 'lego96.py')
		self.lego96Data = os.path.join(self.outDir, 'lego96.dat')
		code = [self.software.python + ' ' + self.software.lego96,
			'--nameVcf ' + self.nameVcf,
			'--reference ' + self.database.reference,
			'--outdir ' + self.outDir
		]
		return '\\\n\t '.join(code)

	def BayesNMF(self):
		self.software.bayesnmf = os.path.join(sys.path[0], 'BayesNMF.signature_discovery.R')
		self.activityBarplot = None
		self.signaturePlot = None
		code=[self.software.R_2_15_3, self.software.bayesnmf, self.lego96Data, self.tumorType, self.outDir]
		return ' '.join(code) + '\n'

	def sigCluster(self):
		pass

		
	def printJob(self):
		code = [self.LeGo96(), self.BayesNMF()]
		return ' &&\n'.join(code)

	def __call__(self):
		subprocess.check_call(self.printJob(), shell=True)
		


def paramsParse():
	parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--nameVcf', help = "Name\tsSNVVCF\nName\tsSNVVCF\n...",required = True)
	parser.add_argument('--outdir', help = "output dir", default = '.')
	parser.add_argument('--tumorType', help = 'cohort name', default = 'Cohort')
	return parser.parse_args()

if __name__ == '__main__':
	args=paramsParse()
	bayesNMF=BAYESNMF(args.nameVcf, args.outdir, args.tumorType)
	bayesNMF()