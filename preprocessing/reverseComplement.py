__author__ = 'pedro'

import argparse
import sys
import subprocess

#This script performs the reverse complement of a list of files and automatically creates the new file in the same directory of the original one.
#It uses seqtk tool to do that.

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def revComp(inputFile):
    files_processed = 0
    seqtk_path = "/opt/tools/seqtk/seqtk seq -r"
    with open(inputFile[0]) as file:

        for line in file:
            if '.fq' in line:
                output_str = line[:-4] + "_revComp.fq"
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

