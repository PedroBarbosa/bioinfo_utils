__author__ = 'pedro'

import argparse
import sys
import subprocess
import os

#This script performs the reverse complement of a list of files and automatically creates the new file in the 'revcomp' directory automatically generated.
#It uses seqtk tool to do that.

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def revComp(inputFile):
    files_processed = 0
    seqtk_path = "/opt/tools/seqtk/seqtk seq -r"
    output_str=""
    with open(inputFile[0]) as file:
	
	if not os.path.exists("revcomp"):
    		os.makedirs("revcomp")		
        for line in file:
	    output_str=""
            if '.fq' in line:
                #reformat file name
		split_list = line.split("/")
		filename=split_list[-1]
		split_list.pop()
		split_list.append("revcomp")
		for i in split_list:
			output_str=output_str + i + "/"
		output_str+=filename
		output_str = output_str[:-4] + "_revComp.fq"
               	
		#run seqtk
		out_file = open(output_str,"w")
                cmd = seqtk_path + " " + line
                ps = subprocess.Popen(cmd, shell=True, stdout=out_file)
                ps.wait()
                out_file.close()
                files_processed += 1
                print "Number of files processed:\t", files_processed
    file.close()




parser = MyParser()

parser.add_argument(dest='inputFile', metavar='input_files', nargs=1,
                    help='File containing the list of files to process.')
args = parser.parse_args()


if __name__ == "__main__":
    revComp(args.inputFile)

