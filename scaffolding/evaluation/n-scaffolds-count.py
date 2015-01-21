__author__ = 'pedro'

import scaffolding.evaluation.extractContigs
import argparse
import sys

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)









if __name__ == "__main__":

        parser = MyParser()
        parser.add_argument(dest='inputFile', metavar='input_files', nargs=1,
                    help='Contigs/Scaffolds fasta file to be processed.')
        parser.add_argument('-n', '--minlength', metavar='', type=int, required=True, help='Minimum length of contigs/scaffolds to process')
        parser.add_argument('-o', '--outputFile', metavar='', required = True, help='Output file')
        args = parser.parse_args()

