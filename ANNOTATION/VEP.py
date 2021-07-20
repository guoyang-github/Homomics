# coding=utf-8
# guoyang
# 封装脚本
import subprocess
import argparse
import os.path
import sys
import os


usage='''
重点是对VEP的安装及使用
perl INSTALL.pl -c ./cache
1. 安装API
2. 安装vep cache文件   注释
3. 安装reference     这一步的fa的bgzip压缩和做index可能会有问题，在初始安装的时候要用好 -a al 和 --NO_HTSLIB （The script will also attempt to install a Perl::XS module, Bio::DB::HTS, for rapid access to bgzipped FASTA files. If this fails, you may add the --NO_HTSLIB flag when running the installer; VEP will fall back to using Bio::DB::Fasta for this functionality）
4. 安装plugin


VEP
The VEP uses "cache files" or a remote database to read genomic data. Using cache files gives the best performance.
By default the VEP installs cache files in a folder in your home area ($HOME/.vep); you can easily change this using the -d flag when running the install script.

The VEP can use FASTA files to retrieve sequence data for HGVS notations and reference sequence checks.

The FASTA file should be automatically detected by the VEP when using --cache or --offline. If it is not, use "--fasta ./cache/homo_sapiens/86_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

Download and install
http://asia.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer

Basic tutorial
http://asia.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html

Running the VEP 参数设置
http://asia.ensembl.org/info/docs/tools/vep/script/vep_options.html

Calculated variant consequences
http://asia.ensembl.org/info/genome/variation/predicted_data.html#consequences
For each variant that is mapped to the reference genome, we identify each Ensembl transcript that overlap the variant. We then use a rule-based approach to predict the effects that each allele of the variant may have on the transcript.

Ensembl Variation - Predicted data
http://asia.ensembl.org/info/genome/variation/predicted_data.html

Ensembl Variation - Data description
http://asia.ensembl.org/info/genome/variation/data_description.html


使用：
variant_effect_predictor.pl
perl variant_effect_predictor.pl -i example_GRCh38.vcf --cache
过滤：
perl filter_vep.pl -i variant_effect_output.txt -filter "SIFT is deleterious" | head -n15
'''


class innerSoftware:
    def __init__(self):
        self.perl='perl'
        self.vep='variant_effect_predictor.pl'

class innerDatabase:
    def __init__(self):
        self.vepData = '/ensembl-tools-release-86/scripts/variant_effect_predictor/cache'
        self.vepReference = '/ensembl-tools-release-86/scripts/variant_effect_predictor/cache/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
	

class VEP:
    def __init__(self, infile, outfile, fork, software = innerSoftware(), database = innerDatabase()):
        self.infile = os.path.abspath(infile)
        self.outfile = os.path.abspath(outfile)
        self.fork = fork
        self.software = software
        self.database = database
		
    def variantEffectPredictor(self):
        code=[self.software.perl + ' ' + self.software.vep,
            ###Basic options###
            #'--everything',   #Shortcut flag to switch on all of the following:   --sift b, --polyphen b, --ccds, --uniprot, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, --tsl, --appris, --gene_phenotype --gmaf, --maf_1kg, --maf_esp, --maf_exac, --pubmed, --variant_class
            '--fork ' + str(self.fork),   #Enable forking, using the specified number of forks. Forking can dramatically improve the runtime of the script. Not used by default

            ###Input options###
            '--species homo_sapiens',   #Species for your data. This can be the latin name e.g. "homo_sapiens" or any Ensembl alias e.g. "mouse". Specifying the latin name can speed up initial database connection as the registry does not have to load all available database aliases on the server. Default = "homo_sapiens"
            '--assembly GRCh37',    #Select the assembly version to use if more than one available. If using the cache, you must have the appropriate assembly's cache file installed. If not specified and you have only 1 assembly version installed, this will be chosen by default. Default = use found assembly version
            '--input_file ' +  self.infile,  #Input file name. If not specified, the script will attempt to read from STDIN.
            '--format vcf',   #Input file format - one of "ensembl", "vcf", "pileup", "hgvs", "id". By default, the script auto-detects the input file format. Using this option you can force the script to read the input file as Ensembl, VCF, pileup or HGVS format, a list of variant identifiers (e.g. rsIDs from dbSNP), or the output from VEP (e.g. to add custom annotation to an existing results file using --custom). Auto-detects format by default
            '--output_file ' + self.outfile,   #Output file name. The script can write to STDOUT by specifying STDOUT as the output file name - this will force quiet mode. Default = "variant_effect_output.txt"
            '--force_overwrite',   #By default, the script will fail with an error if the output file already exists. You can force the overwrite of the existing file by using this flag. Not used by default
            '--html',   #Generate an additional HTML version of the output file containing hyperlinks to Ensembl and other resources. File name of this file is [output_file].html

            ###Cache options###
            '--cache',   #Enables use of the cache. Add --refseq to use the refseq cache (if installed). 
            '--dir ' + self.database.vepData,   #Specify the base cache/plugin directory to use. Default = "$HOME/.vep/"
            '--offline',   #Enable offline mode. No database connections will be made, and only a complete cache (either downloaded or built using --build) can be used for this mode. Add --refseq to use the refseq cache (if installed). Not used by default
            '--fasta ' + self.database.vepReference,   #Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence. The first time you run the script with this parameter an index will be built which can take a few minutes. This is required if fetching HGVS annotations (--hgvs) or checking reference sequences (--check_ref) in offline mode (--offline), and optional with some performance increase in cache mode (--cache). See documentation for more details. Not used by default

            ###Output options###选择注释内容
            '--variant_class',   #Output the Sequence Ontology variant class. Not used by default
            '--sift b',   #p|s|b  prediction term, score or both
            '--polyphen b',
            '--humdiv',  #Human only Retrieve the humDiv PolyPhen prediction instead of the defaulat humVar. Not used by default
            '--gene_phenotype',   #Indicates if the overlapped gene is associated with a phenotype, disease or trait. See list of phenotype sources. Not used by default
            '--regulatory',   #Look for overlaps with regulatory regions. The script can also call if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature. Not used by default
            #'--cell_type',   #Report only regulatory regions that are found in the given cell type(s). Can be a single cell type or a comma-separated list. The functional type in each cell type is reported under CELL_TYPE in the output. To retrieve a list of cell types, use --cell_type list. Not used by default
            #'--custom [filename]'   #Add custom annotation to the output. Files must be tabix indexed or in the bigWig format. Multiple files can be specified by supplying the --custom flag multiple times. See here for full details. Not used by default
            #'--plugin [plugin name]'   #Use named plugin. Plugin modules should be installed in the Plugins subdirectory of the VEP cache directory (defaults to $HOME/.vep/). Multiple plugins can be used by supplying the --plugin flag multiple times. See plugin documentation. Not used by default
            #'--individual [all|ind list]'   #Consider only alternate alleles present in the genotypes of the specified individual(s). May be a single individual, a comma-separated list or "all" to assess all individuals separately. Individual variant combinations homozygous for the given reference allele will not be reported. Each individual and variant combination is given on a separate line of output. Only works with VCF files containing individual genotype data; individual IDs are taken from column headers. Not used by default
            '--allele_number',   #Identify allele number from VCF input, where 1 = first ALT allele, 2 = second ALT allele etc. Not used by default
            '--total_length',   #Give cDNA, CDS and protein positions as Position/Length. Not used by default
            '--numbers',   #Adds affected exon and intron numbering to to output. Format is Number/Total. Not used by default
            '--domains',   #Adds names of overlapping protein domains to output. Not used by default
            '--no_escape',   #Don't URI escape HGVS strings. Default = escape
            #'--terms [ensembl|so]'   #The type of consequence terms to output. The Ensembl terms are described here. The Sequence Ontology is a joint effort by genome annotation centres to standardise descriptions of biological sequences. Default = "SO"

            ###Identifiers###
            '--hgvs',   #Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate. To generate HGVS identifiers when using --cache or --offline you must use a FASTA file and --fasta. HGVS notations given on Ensembl identifiers are versioned. Not used by default
            '--shift_hgvs 1',   #Enable or disable 3' shifting of HGVS notations. When enabled, this causes ambiguous insertions or deletions (typically in repetetive sequence tracts) to be "shifted" to their most 3' possible coordinates (relative to the transcript sequence and strand) before the HGVS notations are calculated; the flag HGVS_OFFSET is set to the number of bases by which the variant has shifted, relative to the input genomic coordinates. Disabling retains the original input coordinates of the variant. Default: 1 (shift)
            '--protein',   #Add the Ensembl protein identifier to the output where appropriate. Not used by default
            '--symbol',   #Adds the gene symbol (e.g. HGNC) (where available) to the output. Not used by default
            '--ccds',   #Adds the CCDS transcript identifer (where available) to the output. Not used by default
            '--uniprot',   #Adds best match accessions for translated protein products from three UniProt-related databases (SWISSPROT, TREMBL and UniParc) to the output. Not used by default
            '--tsl',   #Adds the transcript support level for this transcript to the output. Not used by default
            '--appris',   #Adds the APPRIS isoform annotation for this transcript to the output. Not used by default
            '--canonical',   #Adds a flag indicating if the transcript is the canonical transcript for the gene. Not used by default
            '--biotype',   #Adds the biotype of the transcript or regulatory feature. Not used by default
            '--xref_refseq',   #Output aligned RefSeq mRNA identifier for transcript. NB: theRefSeq and Ensembl transcripts aligned in this way MAY NOT, AND FREQUENTLY WILL NOT, match exactly in sequence, exon structure and protein product. Not used by default

            ###Co-located variants###
            '--check_existing',   #Checks for the existence of known variants that are co-located with your input. By default the alleles are not compared - to do so, use --check_alleles. Not used by default
            '--check_alleles',   #When checking for existing variants, only report a co-located variant if none of the alleles supplied are novel. For example, if the user input has alleles A/G, and an existing co-located variant has alleles A/C, the co-located variant will not be reported. Strand is also taken into account - in the same example, if the user input has alleles T/G but on the negative strand, then the co-located variant will be reported since its alleles match the reverse complement of user input. Not used by default
            #'--check_svs',   #Checks for the existence of structural variants that overlap your input. Currently requires database access. Not used by default
            '--gmaf',   #Add the global minor allele frequency (MAF) from 1000 Genomes Phase 3 data for any existing variant to the output. Not used by default
            '--maf_1kg',   #Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output. Note the reported allele(s) and frequencies are for the non-reference allele from the original data, not necessarily the alternate allele from user input. Must be used with --cache Not used by default
            '--maf_esp',   #Include allele frequency from NHLBI-ESP populations. Note the reported allele(s) and frequencies are for the non-reference allele from the originial data, not necessarily the alternate allele from user input. Must be used with --cache Not used by default
            '--maf_exac',   #Include allele frequency from ExAC project populations. Note the reported allele(s) and frequencies are for the non-reference allele from the originial data, not necessarily the alternate allele from user input. Must be used with --cache Not used by default
            '--pubmed',   #Report Pubmed IDs for publications that cite existing variant. Must be used with --cache. Not used by default
            '--failed 1',   #When checking for co-located variants, by default the script will exclude variants that have been flagged as failed. Set this flag to include such variants. Default: 0 (exclude)

            ###Data format options###
            #'--vcf',   #Writes output in VCF format. Consequences are added in the INFO field of the VCF file, using the key "CSQ". Data fields are encoded separated by "|"; the order of fields is written in the VCF header. Output fields can be selected by using --fields. If the input format was VCF, the file will remain unchanged save for the addition of the CSQ field (unless using any filtering). Custom data added with --custom are added as separate fields, using the key specified for each data file. Commas in fields are replaced with ampersands (&) to preserve VCF format. Not used by default
            '--tab',   #Writes output in tab-delimited format. Not used by default
            #'--json',   #Writes output in JSON format. Not used by default
            #'--gvf',   #Writes output in GVF format. Not used by default
            #'--fields [list]'   #Configure the output format using a comma separated list of fields. Fields may be those present in the default output columns, or any of those that appear in the Extra column (including those added by plugins or custom annotations). Output remains tab-delimited. Not used by default
            '--minimal',   #Convert alleles to their most minimal representation before consequence calculation i.e. sequence that is identical between each pair of reference and alternate alleles is trimmed off from both ends, with coordinates adjusted accordingly. Note this may lead to discrepancies between input coordinates and coordinates reported by VEP relative to transcript sequences; to avoid issues, use --allele_number and/or ensure that your input variants have unique identifiers. The MINIMISED flag is set in the VEP output where relevant. Not used by default

            ###Filtering and QC options###
            '--check_ref',   #Force the script to check the supplied reference allele against the sequence stored in the Ensembl Core database. Lines that do not match are skipped. Not used by default
            #'--pick',   #Pick once line or block of consequence data per variant, including transcript-specific columns. Consequences are chosen according to the criteria described here, and the order the criteria are applied may be customised with --pick_order. This is the best method to use if you are interested only in one consequence per variant. Not used by default
            #this is the option we anticipate will be of use to most users. VEP chooses one block of annotation per variant, using an ordered set of criteria. This order may be customised using --pick_order.
            #    1.canonical status of transcript
            #    2.APPRIS isoform annotation
            #    3.transcript support level
            #    4.biotype of transcript (protein_coding preferred)
            #    5.CCDS status of transcript
            #    6.consequence rank according to this table
            #    7.translated, transcript or feature length (longer preferred)
            #Note that some categories may not be available for the species or cache version that you are using; in these cases the category will be skipped and the next in line used. 
            #'--flag_pick',   #As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
            #'--pick_allele',   #Like --pick, but chooses one line or block of consequence data per variant allele. Will only differ in behaviour from --pick when the input variant has multiple alternate alleles. Not used by default
            '--flag_pick_allele',   #As per --pick_allele, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
            '--pick_order canonical,appris,tsl,biotype,ccds,rank,length'   #Customise the order of criteria applied when choosing a block of annotation data with e.g. --pick. See this page for the default order. Valid criteria are: canonical,appris,tsl,biotype,ccds,rank,length
        ]		
        return '   \\\n\t'.join(code)
		
		
    def printJob(self):
        return self.variantEffectPredictor() + '\n\n'

    def __call__(self):
        subprocess.check_call(self.printJob(), shell=True)
		


def paramsParse():
    parser=argparse.ArgumentParser(description=usage, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--infile',help = "input vcf file",required = True)
    parser.add_argument('--outfile',help = "output vep-annoted vcf／tab file",required = True)
    parser.add_argument('--fork',help = "threads to use", default = 4, type = int)

    return parser.parse_args()

if __name__ == '__main__':
    args = paramsParse()
    vep = VEP(args.infile, args.outfile, args.fork)
    vep()