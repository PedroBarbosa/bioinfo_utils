from Bio import SeqIO
import argparse
import os 
def transformFastq(inputF):
	for file in inputF:
		print("Processing %s file.." % file)
		sequences = SeqIO.parse(file, 'fastq') 
		#for record in sequences:
		basename = os.path.basename(file).rsplit('.',1)[0]
		SeqIO.convert(file, "fastq", os.getcwd() +  "/" + basename + ".fasta" , "fasta")
		SeqIO.convert(file, "fastq", os.getcwd() +  "/" + basename + ".qual", "qual")
		print("Done")
parser = argparse.ArgumentParser(description = 'Simple script to convert fastq file/s to fasta and qual. Outputs filenames will be based on the input file, except extension.')
parser.add_argument(dest='inputFile', metavar='input_files', nargs= '+', help='Fastq file/s to be processed.')
args = parser.parse_args()

transformFastq(args.inputFile)
