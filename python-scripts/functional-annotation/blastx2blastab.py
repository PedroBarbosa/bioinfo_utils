import argparse
import re

__author__ = 'pedro'

parser = argparse.ArgumentParser(description='This is a script to parse blastx files into blasttab')
parser.add_argument('-i', dest='blastx_results_file', help='Blast results file', required=True)
parser.add_argument('-o', dest='outputed_file', help='Name of the BLAST tabular output file', required=True)
args = parser.parse_args()


def parseblastXfile(blastFile):
    print "Processing file ..."
    filename = args.outputed_file
    fileoutput = open (filename, 'w')
    fields = "#Fields:  Query   Subject identity    aln-len mismatch    gap-openings    q.start q.end   s.start s.end   e-value bit-score\n"
    fileoutput.write(fields)
    new_hit = [None] * 12
    first_hit = True
    with open(blastFile) as file:
        for line in file:

            line = line.strip()
            if line.startswith('Query='):
                query = ""
                query = line.split('Query=',1)[1]  #get query
                first_hit = True

            elif line.startswith('>'):
                new_hit = [None] * 12
                new_hit[0] = query   ##add query
                subject = re.findall(r">(\S+).*", line) # r"# line.split('>',1)[1].strip() #get subject
                new_hit[1] = subject[0]       ##add subject

            elif "Score" in line:
                metrics = re.findall(r"(?:\d+(?:[eE\.][-]?\d*))", line)
                if len(metrics) is 1:
                    new_hit[10] = line[-1:]  ##evalue when is positive (Very bad)
                    new_hit[11] = metrics[0] #score
                else:
                    new_hit[11] = metrics[0] #score
                    new_hit[10] = metrics[1] #evalue


            elif "Identities" in line:

                alignment_features = line.split(',')
                identity = re.findall("\((\d+).*", alignment_features[0]) #get identity
                aln_len = re.findall("\/(\d+).*", alignment_features[0])  #get alignment length
                mismatch = alignment_features[0].split("/")  #get number of mismatches
                mismatch_value = int(mismatch[1].replace(mismatch[1],mismatch[1][:2])) - int(mismatch[0].replace(mismatch[0],mismatch[0][-2:]))
                gaps = re.findall("(\d+)\/.*",alignment_features[2]) #get number of gaps
                new_hit[2] = identity[0]
                new_hit[3] = aln_len[0]
                new_hit[4] = str(mismatch_value)
                new_hit[5] = gaps[0]

            elif line.startswith('Query:'):

                query_coord = re.findall(r'\d+', line)
                new_hit[6] = query_coord[0]    ## add query start
                new_hit[7] = query_coord[1]    ## add query end

            elif line.startswith('Sbjct:'):

                subject_coord = re.findall(r'\d+', line)
                new_hit[8] = subject_coord[0]    ## add subject start
                new_hit[9] = subject_coord[1]    ## add subject end
                fileoutput.write('\t'.join(new_hit) + '\n')

    return fileoutput

parseblastXfile(args.blastx_results_file)