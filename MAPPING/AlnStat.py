# coding=utf-8
#guoyang 2016.5.19
import HTSeq
import argparse
import os.path


usage='''
Flagstat

flags = {0x1 : "template having multiple segments in sequencing",
         0x2 : "each segment properly aligned according to the aligner",
         0x4 : "segment unmapped",
         0x8 : "next segment in the template unmapped",
         0x10 : "SEQ being reverse complemented",
         0x20 : "SEQ of the next segment in the template being reversed",
         0x40 : "the first segment in the template",
         0x80 : "the last segment in the template",
         0x100 : "secondary alignment",
         0x200 : "not passing quality controls",
         0x400 : "PCR or optical duplicate",
         0x800 : "supplementary alignment"}
		 


		 
Total:	108447138 (100%)
Duplicate:	18349819 (16.94%)
Mapped:	108309617 (99.87%)
Properly mapped:	106895214 (98.57%)
PE mapped:	108191438 (99.76%)
SE mapped:	236358 (0.22%)
With mate mapped to a different chr:	1096480 (1.01%)
With mate mapped to a different chr ((mapQ>=5)):	965460 (0.89%)




108661519 in total
0 QC failure
18349819 duplicates
108523998 mapped (99.87%)
108661519 paired in sequencing
54323372 read1
54338147 read2
107012199 properly paired (98.48%)
108404106 with itself and mate mapped
119892 singletons (0.11%)
1226788 with mate mapped to a different chr
1001006 with mate mapped to a different chr (mapQ>=5)
'''



def paramsParse():
	parser = argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('--bam',help="the input bam file", required=True)
	parser.add_argument('--outdir',help="output dir", required=True)
	return parser.parse_args()


def flagstat(bamfile):
	mate_mapped_same_chr=0
	mate_mapped_dif_chr=0
	mate_mapped_dif_chr_a5=0
	mapped=0
	unmapped=0
	paired=0
	singletons=0
	read1=0
	read2=0
	properly=0
	duplicate=0
	total=0
	output=[]
	
	bamFile = HTSeq.BAM_Reader(bamfile)
	for almnt in bamFile:
		if almnt.aligned:  #mapped read
			if almnt.not_primary_alignment: #if it is not the primaray alignment   这应该就是和GATK及samtools的flagstat不同的地方，我们只考虑primary比对
				continue
			total += 1
			mapped += 1
			if almnt.mate_aligned: 
				paired +=1
				if almnt.proper_pair:
					properly += 1
				if almnt.iv.chrom == almnt.mate_start.chrom:
					mate_mapped_same_chr += 1
				else:
					mate_mapped_dif_chr += 1
					if almnt.aQual >= 5:
						mate_mapped_dif_chr_a5 += 1
			else:
				singletons+=1
			if almnt.pe_which == 'first':
				read1 += 1
			if almnt.pe_which == 'second':
				read2 += 1
			if almnt.pcr_or_optical_duplicate:
				duplicate+=1
		else:
			unmapped += 1
			total += 1
	
	#注意递进的比例关系					
	output.append('Total:\t' + str(total) + ' (100%)')
	output.append('Mapped:\t' + str(mapped) + ' (' + str(float(mapped)/total*100) + '%)')
	output.append('Duplicate:\t' + str(duplicate) + ' (' + str(float(duplicate)/mapped*100) + '%)')
	output.append('PE mapped:\t' + str(paired) + ' (' + str(float(paired)/mapped*100) + '%)')
	output.append('Properly mapped:\t' + str(properly) + ' (' + str(float(properly)/paired*100) + '%)')
	output.append('With mate mapped to a different chr:\t' + str(mate_mapped_dif_chr) + ' (' + str(float(mate_mapped_dif_chr)/paired*100) + '%)')
	output.append('With mate mapped to a different chr ((mapQ>=5)):\t' + str(mate_mapped_dif_chr_a5) + ' (' + str(float(mate_mapped_dif_chr_a5)/paired*100) + '%)')
	output.append('SE mapped:\t' + str(singletons) + str(float(singletons)/mapped*100) + '%)')
	output.append('Read1:\t' + str(read1))
	output.append('Read2:\t' + str(read2))
	
	return '\n'.join(output)
	
	
	
if __name__ == '__main__':
	args=paramsParse()
	outfile=os.path.join(os.path.abspath(args.outdir), os.path.basename(args.bam).replace('.bam','') + '.flagstat')
	with open(outfile,'w') as output:
		output.write(flagstat(args.bam))

