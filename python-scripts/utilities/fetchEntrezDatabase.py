import argparse
import os
parser = argparse.ArgumentParser(description='Script to get a list of accession numbers from a given taxonomic ID and all its sub txids. Need Biopython installed, a python27 '
                                             ' environment and a file mapping the old GIs to new accesssion numbers (downloadable through the NCBI ftp site). Also requires the'
                                             ' script gi2acession.py script from the same NCBI directory.')
parser.add_argument(dest='output', metavar='outputBasename', help='Output basename to write the results.')
parser.add_argument('-e', metavar='email', required="True", help='Provide an email to access the NCBI Entrez database.')
parser.add_argument('-db', metavar='database', default='protein', help='Database to search the given terms. Default:"protein".')
parser.add_argument('-t',metavar='taxid', default='33090',help='Top level Taxonomic ID to search all accession number. Default:"Green plants: 33090"')
parser.add_argument('-p', metavar='pathScript',default=os.getcwd()+"/gi2accession.py", help='Path to the auxiliar script to convert gi2accession. Default:this script directory.')
parser.add_argument('-d', '--download', action='store_true', help='If no auxiliar files are available, automatic download all the data.')
args=parser.parse_args()

#import importlib.util
#spec = importlib.util.spec_from_file_location("gi2accession", args.p)
#exec = importlib.util.module_from_spec(spec)
from Bio import Entrez

def getGIsFromTxid(email,taxid,db,outbasename):
    Entrez.email = email
    #record = Entrez.read(Entrez.einfo())
    print("Retrieving list of old GIs from Entrez database..")
    handle = Entrez.esearch(db=db, term="txid"+taxid+"[Organism]",retmax=10**9)
    list=Entrez.read(handle)['IdList']
    out_gi = outbasename + '_gi.txt'
    with open(out_gi, 'w') as outfile:
        outfile.write('\n'.join(list))
    return out_gi


def getAccFromGIfile(out_gi,outbasename,path):
    out_acc = outbasename + "_acession.txt"
    #print(exec)
    #spec.loader.exec_module(exec)
    print("Extracting accession numbers..")
    os.system("python " + path + " < " + out_gi + " > " + out_acc)
    print("Done.")
if not args.download:
    str_out =getGIsFromTxid(args.e,args.t,args.db,args.output)
    getAccFromGIfile(str_out, args.output,args.p)
else:
    print("Function not ready. Download by yourself.")