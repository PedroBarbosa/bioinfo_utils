import argparse
import sys
import os
import subprocess
from collections import defaultdict

__author__ = 'pedro'


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        self.print_help()
        sys.exit(2)

#parser = argparse.ArgumentParser(description='This is a script to compare blast tabular results with different search engines.')
parser = MyParser()
parser.add_argument(dest='inputFiles', metavar='annotated_files', nargs='+',
                    help='Annotation files to be analyzed (minimum 2).')
parser.add_argument('-r', '--reads', type=argparse.FileType('r'), metavar='', help='Original reads file')
parser.add_argument('-s', '--subsample', action='store_true',
                    help='If true, the reads file in FASTA should be provided with the -r argument'
                         ' and the number of reads to be sampled passed with the -n argument. Default: "False".')
parser.add_argument('-n', '--subsample_size', type=int, metavar='',
                    help='Number of reads to be sampled from the original reads file. Use only if" -s" is True.')
parser.add_argument('-k', '--kegg', action='store_true',
                    help='Flag to indicate if the annotation was performed against kegg. If True, some additional'
                         ' statistics will be generated.')
parser.add_argument('-b', '--blast', action='store_true',
                    help='Is one of the files a true blast annotation? If so, place it as the first positional '
                         'argument.')
args = parser.parse_args()


def createDict():
    dict = defaultdict(list)
    return dict



def processAnnotatedReadsBlastTab(blastFiles, dict, numb_files):
    print 'Processing files..'

    for i in range(0, num_files):

        #process = subprocess.Popen(['wc','-l'], stdin=blastFiles[i], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        n = int(subprocess.check_output(['wc', '-l', blastFiles[i]]).split()[0]) + 1   #number of lines

        previous_line = ""
        numb_hits = -1
        iter = 0
        if i == 1:
            print 'Finished processing first file..'
        with open(blastFiles[i]) as file:
            for line in file:

                iter += 1
                if not line.startswith('#'):
                    read = line.split("\t")[0]  #read itself

                    if read == previous_line:
                        numb_hits += 1

                    else:
                        if numb_hits != -1:
                            update_list = dict.get(previous_line)
                            update_list[i * 3 + 1] = numb_hits + 1  #position of the nnumber of hits.
                            dict[previous_line] = update_list

                        listOfValues = [None] * 3 * num_files
                        numb_hits = 0
                        previous_line = read

                        if read not in dict:
                            listOfValues[i * 3 + 0] = 1
                            dict[read] = listOfValues


                        #elif dict[read].count(None) != len(listOfValues):  ## update hash table for files other than the first
                        else:
                            temp = dict.get(read)
                            temp[i * 3 + 0] = 1
                            dict[read] = temp


                    if iter == n:
                        update_list = dict.get(previous_line)
                        update_list[i * 3 + 1] = numb_hits + 1  #position of the number of hits.
                        dict[previous_line] = update_list

    return dict


def getStats2samples(processed_dict):
    print "Calculating stats .."
    both_present = only_first = only_second = 0
    unique_reads_first = unique_reads_second = []
    both_significant = both_one = significant_first = significant_second = one_hit_first = one_hit_second = 0
    for key, value in processed_dict.iteritems():

        if value[0] and value[3] == 1:
            both_present += 1
        elif value[0] == 1 and value[3] == None:
            only_first +=1
            unique_reads_first.append(key)
        else:
            only_second +=1
            unique_reads_second.append(key)


        if value[1] and value[4] >= 20:
            both_significant+=1

        elif  value[1] and value[4] == 1:
            both_one +=1

        elif value[1] >= 20:
            significant_first +=1

        elif value[4] >= 20:
            significant_second+=1


    print "Number of reads annotated in both annotations:\t" , both_present
    print "Number of reads annotated only in the first file:\t" ,  only_first
    print "Number of reads annotated only in the second file:\t" , only_second, "\n"
    print "Number of reads annotated with significant hits in both annotations:\t" , both_significant
    print "Number of reads annotated with significant hits only in first file:\t" ,  significant_first
    print "Number of reads annotated with significant hits only in the second file:\t" , significant_second
    print "Number of reads annotated with only one hit in both annotations:\t" , both_one
    samplename = os.getcwd().split("/")[-1]
    filename = os.getcwd() + "/" + samplename + "_blastcomparison.txt"
    fileoutput = open (filename, 'w')
    fileoutput.write("Number of reads annotated in both annotations:\t" + str(both_present) + "\n")
    fileoutput.write("Number of reads annotated only in the first file:\t" + str(only_first) + "\n")
    fileoutput.write("Number of reads annotated only in the second file:\t" + str(only_second) + "\n")
    fileoutput.write("Number of reads annotated with significant hits in both annotations:\t" + str(both_significant) + "\n")
    fileoutput.write("Number of reads annotated with significant hits only in first file:\t" +  str(significant_first) + "\n")
    fileoutput.write("Number of reads annotated with significant hits only in the second file:\t" + str(significant_second) + "\n")
    fileoutput.write("Number of reads annotated with only one hit in both annotations:\t" + str(both_one) + "\n")

    fileoutput.close()


###############################################################################
########################## COMMAND LINE PARSING ###############################
##############################################################################
if len(args.inputFiles) < 2 :  ##test for file
    sys.stderr.write('Error: %s\n' % 'You should specify at least two annotation files.')
    sys.exit(2)

for file in args.inputFiles:
    if not os.path.exists(file):
        print 'Error: Can not open', file, 'file'
        sys.exit((2))

if not args.reads:
    sys.stderr.write('Error: %s\n' % 'Please provide -r argument.')
    sys.exit((2))

if args.subsample:
    if not args.subsample_size:
        sys.stderr.write('Error: %s\n' % 'Please provide -n argument.')
        sys.exit((2))
    elif args.blast:
        sys.stderr.write('Error: %s\n' % 'Please remove -b argument. Not compatible with -s.')
        sys.exit(2)
    else:
        print 'ya'
        #Do stuff
elif args.subsample_size:
    sys.stderr.write('Error: %s\n' % 'Please remove -n argument. Only required when -s is provided.')
    sys.exit(2)

elif args.blast:

    for firstline in args.inputFiles[0]:
        if firstline.startswith('blastx'):
            print 'Blast file ok.'
            break
        else:
            sys.stderr.write(
                'Error: %s\n' % 'It seems that the first annotated file provided does not come from a ncbi blast run. Please provide '
                                'a blast file as the first file.')
            sys.exit(2)

            #do stuff

else:
    num_files = len(args.inputFiles)
    hash = createDict()
    processed_hash = processAnnotatedReadsBlastTab(args.inputFiles, hash, num_files)
    getStats2samples(processed_hash)

    #for key, list in processed_hash.iteritems():
     #   print key, list
        #for value in list:
        #    print('{} : {}'.format(value))
##useful comparison from rapsearch: The overlap lists the number of reads that have homologs in several search engines. Only one of them.
#Also useful to allow subsample some blast results of the reads that have hits only in one search engine

#use measure to compare reads with only one hit and reads with significant hits >20 ??

# another measure of hits assignments - the percentage of hits that account for the same KO (when searching agains kegg)


