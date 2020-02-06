from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Script to creAte fasta aminoacid file from genbank using CDS seq feature types..')
parser.add_argument(dest='genbank', metavar='genbank', help='Genbank file to process.')
parser.add_argument(dest='fasta', metavar='fasta', help='Fasta output file.')
parser.add_argument('-n','--noTranlation', action="store_true", help='If genbank record doesn\'t have translation qualifier, set this option to do it on the fly.')
args=parser.parse_args()


gbk_filename = args.genbank
faa_filename = args.fasta
input_handle  = open(gbk_filename, "r")
output_handle = open(faa_filename, "w")
g_id = 1
for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" and not args.noTranlation:
            try:
                assert len(seq_feature.qualifiers['translation'])==1
                output_handle.write(">%s from %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   seq_feature.qualifiers['translation'][0]))
            except KeyError:
                print("No translation qualifier found in genbank record. Please set --noTranslation flag on the command line.")
                exit(1)


        elif seq_feature.type=="CDS" and args.noTranlation:

            if seq_feature.location.strand == 1:
                pep = seq_record.seq[seq_feature.location.start:seq_feature.location.end].translate(table=11)
                #dna2 = seq_record.seq[seq_feature.location.start:seq_feature.location.end].translate(table=11,to_stop=True) #truncate when finds stop codon
                try:
                    output_handle.write(">MAPK_%i %s %s\n" % (g_id,seq_feature.qualifiers['protein_id'][0],seq_feature.qualifiers['product'][0]))
                    output_handle.write("%s\n" % str(pep))
                except:
                    output_handle.write(">MAPK_%i\n" % g_id)
                    output_handle.write("%s\n" % str(pep))
                #print("feature\t" + str(len(seq_feature)))

            if seq_feature.location.strand == -1:
                rev = seq_record.seq[seq_feature.location.start:seq_feature.location.end].reverse_complement().translate(table=11)
                output_handle.write(">MAPK_%i %s %s\n" % (g_id,seq_feature.qualifiers['protein_id'][0],seq_feature.qualifiers['product'][0]))
                output_handle.write("%s\n" % str(rev))
                #print(seq_feature.location.start)
            g_id += 1
output_handle.close()
input_handle.close()
print("Done")