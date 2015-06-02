__author__ = 'pedro'

import argparse
import os
import csv
import numpy as np
from collections import defaultdict

#0 - Protein Accession (e.g. P51587)
#1 - Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
#2   Sequence Length (e.g. 3418)
#3   Analysis (e.g. Pfam / PRINTS / Gene3D)
#4   Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
#5   Signature Description (e.g. BRCA2 repeat profile)
#6   Start location
#7   Stop location
#8   Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
#9   Status - is the status of the match (T: true)
#10  Date - is the date of the run
#11  (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprscan option is switched on)
#12  (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprscan option is switched on)
#13  (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
#14  (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)


#AVAILABLE ANALYSIS
#ProDom (2006.1) : ProDom is a comprehensive set of protein domain families automatically generated from the UniProt Knowledge Database.
#Hamap (201502.04) : High-quality Automated and Manual Annotation of Microbial Proteomes
#SMART (6.2) : SMART allows the identification and analysis of domain architectures based on Hidden Markov Models or HMMs
#ProSiteProfiles (20.105) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
#ProSitePatterns (20.105) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them
#SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotation for all proteins and genomes.
#PRINTS (42.0) : A fingerprint is a group of conserved motifs used to characterise a protein family
#PANTHER (9.0) : The PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System is a unique resource that classifies genes by their functions, using published scientific experimental evidence and evolutionary relationships to predict function even in the absence of direct experimental evidence.
#Gene3D (3.5.0) : Structural assignment for whole genes and genomes using the CATH domain structure database
#PIRSF (3.01) : The PIRSF concept is being used as a guiding principle to provide comprehensive and non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect their evolutionary relationships.
#Pfam (27.0) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs)
#TIGRFAM (15.0) : TIGRFAMs are protein families based on Hidden Markov Models or HMMs
#Coils (2.2)


def processInterproRun(interproScanFile,):
        hashmap = defaultdict(list)
        with open(interproScanFile) as file:
            for line in file:
                attributes = line.rstrip().split("\t")
                # attributes=filter(None,attributes)

                if not line.startswith('#'):
                    gene = attributes[0]
                    base_tuple = (attributes[3], attributes[4], attributes[5], attributes[8])

                    if len(attributes) > 11 and len(attributes) < 14:  #with InterPro annotations
                        tuple_interpro = base_tuple + (attributes[11],) + (attributes[12],)
                        hashmap[gene].append(tuple_interpro)

                    elif len(attributes) == 14:  #with GO annotations
                        tuple_go = base_tuple + (attributes[11],) + (attributes[12],) + (attributes[13],)
                        hashmap[gene].append(tuple_go)

                    elif len(attributes) == 15:  #with pathway information
                        tuple_go_path = base_tuple + (attributes[11],) + (attributes[12],) +(attributes[13],) + (attributes[14],)
                        hashmap[gene].append(tuple_go_path)
                    else:
                        hashmap[gene].append(base_tuple)

        file.close()
        return hashmap



class InterproDomains:
    def __init__(self):
       self.genes_with_domains=self.total_number_domains= self.prodom_hits= self.hamap_hits= self.smart_hits= self.prositeProfiles_hits= self.prositePatterns_hits= \
       self.superfamily_hits= self.prints_hits= self.panther_hits= self.gene3D_hits= self.pirsf_hits= self.pfam_hits= self.tigrfam_hits= self.coils_hits=\
       self.interpro_hits=self.go_terms=self.pathway_annotation=self.kegg_pathways= self.kegg_enzymes = self.metacyc = self.reactome = self.unipathway  = 0

       self.unique_prodom_hits= set()
       self.unique_hamap_hits= set()
       self.unique_smart_hits= set()
       self.unique_prositeProfiles_hits = set()
       self.unique_prositePatterns_hits = set()
       self.unique_superfamily_hits = set()
       self.unique_prints_hits = set()
       self.unique_panther_hits = set()
       self.unique_gene3D_hits = set()
       self.unique_pirsf_hits = set()
       self.unique_pfam_hits = set()
       self.unique_tigrfam_hits = set()
       self.unique_interpro_id = set()
       self.unique_go_terms = set()

       self.unique_kegg_pathways = set()
       self.unique_kegg_enzymes = set()
       self.unique_metacyc_pathways = set()
       self.unique_reactome = set ()
       self.unique_unipathway = set()

       self.average_domains_gene = []

    def analyseHash(self,inputDict):

        self.genes_with_domains = len(inputDict)
        for gene,hits in inputDict.iteritems():
            self.average_domains_gene.append(len(hits)) #number of hits for this gene

            #each hit has the following fields: analysis[0],signature_id[1], signature description[2], score[3], interpro_id[4], interpro_description[5], go_terms[6]
            #pathway_analysis[7]
            for hit in hits:
                #print (hit)
                self.total_number_domains +=1


                if len(hit) > 4 and len(hit) < 7: #with interpro annotation
                    self.interpro_hits +=1
                    self.unique_interpro_id.add(hit[4]) #interpro id

                elif len(hit) == 7: #with go annotation
                    self.interpro_hits +=1
                    self.unique_interpro_id.add(hit[4]) #interpro id
                    self.go_terms += 1
                    self.unique_go_terms.add(hit[6]) #go id

                elif len(hit) > 7 : #with pathway annotation
                    self.interpro_hits +=1
                    self.unique_interpro_id.add(hit[4]) #interpro id
                    self.go_terms += 1
                    self.unique_go_terms.add(hit[6]) #go id
                    self.pathway_annotation +=1

                    ##parse pathway field##
                    #print (hit[7])
                    diff_info = hit[7].split("|")
                    for analysis in diff_info:

                        if "KEGG:" in analysis:
                            #pathway
                            kegg=analysis[5:].strip()
                            self.kegg_pathways += 1
                            self.unique_kegg_pathways.add(kegg.split('+')[0])
                            #enzymes
                            for enzyme in kegg.split('+')[1:]:
                                self.kegg_enzymes += 1
                                self.unique_kegg_enzymes.add(enzyme)

                        elif "MetaCyc:" in analysis:
                            metacyc_id=analysis[8:].strip()
                            self.metacyc += 1
                            self.unique_metacyc_pathways.add(metacyc_id)

                        elif "Reactome:" in analysis:
                            reactome_id=analysis[9:].strip()
                            self.reactome += 1
                            self.unique_reactome.add(reactome_id)

                        elif "UniPathway:" in analysis:
                            unipathway_id = analysis[11:].strip()
                            self.unipathway +=1
                            self.unique_unipathway.add(unipathway_id)



                if hit[0] == 'ProDom': #Analysis
                    self.prodom_hits +=1
                    self.unique_prodom_hits.add(hit[1]) #Signature acession

                elif hit[0] == 'Hamap':
                    self.hamap_hits +=1
                    self.unique_hamap_hits.add(hit[1])

                elif hit[0] == 'SMART':
                    self.smart_hits +=1
                    self.unique_smart_hits.add(hit[1])

                elif hit[0] == 'ProSiteProfiles':
                    self.prositeProfiles_hits +=1
                    self.unique_prositeProfiles_hits.add(hit[1])

                elif hit[0] == 'ProSitePatterns':
                    self.prositePatterns_hits +=1
                    self.unique_prositePatterns_hits.add(hit[1])

                elif hit[0] == 'SUPERFAMILY':
                    self.superfamily_hits +=1
                    self.unique_superfamily_hits.add(hit[1])

                elif hit[0] == 'PRINTS':
                    self.prints_hits +=1
                    self.unique_prints_hits.add(hit[1])

                elif hit[0] == 'PANTHER':
                    self.panther_hits +=1
                    self.unique_panther_hits.add(hit[1])

                elif hit[0] == 'Gene3D':
                    self.gene3D_hits +=1
                    self.unique_gene3D_hits.add(hit[1])

                elif hit[0] == 'PIRSF':
                    self.pirsf_hits +=1
                    self.unique_pirsf_hits.add(hit[1])

                elif hit[0] == 'Pfam':
                    self.pfam_hits +=1
                    self.unique_pfam_hits.add(hit[1])

                elif hit[0] == 'TIGRFAM':
                    self.tigrfam_hits +=1
                    self.unique_tigrfam_hits.add(hit[1])

                elif hit[0] == 'Coils':
                    self.coils_hits +=1




    def writeReport(self,outputFile):
            if os.path.exists(outputFile):
                os.remove(outputFile)
            with open(outputFile, "w") as csvfile:
                writer = csv.writer(csvfile,dialect=csv.excel_tab)

                writer.writerow(('Number of genes with a match to an interpro domain:', self.genes_with_domains))
                writer.writerow(('Total number of domains found:', self.total_number_domains))
                writer.writerow(('Average number of matches per gene:', round(np.mean(self.average_domains_gene),4)))
                writer.writerow('')
                writer.writerow(('Total number of domains with Interpro ID:', self.interpro_hits))
                writer.writerow(('Number of different interpro IDs found:', len(self.unique_interpro_id)))
                writer.writerow('')
                writer.writerow(('Total number of domains with GO terms:', self.go_terms))
                writer.writerow(('Number of different GO terms found:', len(self.unique_go_terms)))
                writer.writerow('')
                writer.writerow(('Total number of matches with pathway analysis:', self.pathway_annotation))
                writer.writerow(('','Total number of matches to KEGG pathways:', self.kegg_pathways))
                writer.writerow(('','Total number of matches to KEGG enzymes:', self.kegg_enzymes))
                writer.writerow(('','Total number of matches to MetaCyc database:', self.metacyc))
                writer.writerow(('','Total number of matches to Reactome database:', self.reactome))
                writer.writerow(('','Total number of matches to UniPathway database:', self.unipathway))
                writer.writerow((''))
                writer.writerow(('','Number of different KEGG pathways identified:', len(self.unique_kegg_pathways)))
                writer.writerow(('','Number of different KEGG enzymes identified:', len(self.unique_kegg_enzymes)))
                writer.writerow(('','Number of different MetaCyc matches identified:', len(self.unique_metacyc_pathways)))
                writer.writerow(('','Number of different Reactome mathes identified:', len(self.unique_reactome)))
                writer.writerow(('','Number of different UniPathway matches identified:', len(self.unique_unipathway)))
                writer.writerow('')
                writer.writerow('')
                writer.writerow('')
                writer.writerow('')
                writer.writerow(('Total number of ProDom protein domain families:', self.prodom_hits))
                writer.writerow(('Number of different ProDom protein domains found:', len(self.unique_prodom_hits)))
                writer.writerow('')
                writer.writerow(('Total number of HaMap manually curated family profiles:', self.hamap_hits))
                writer.writerow(('Number of different HaMap manually curated family profiles found:', len(self.unique_hamap_hits)))
                writer.writerow('')
                writer.writerow(('Total number of SMART domain architectures:', self.smart_hits))
                writer.writerow(('Number of different SMART domain architectures found:', len(self.unique_smart_hits)))
                writer.writerow('')
                writer.writerow(('Total number of ProSiteProfiles protein domain families:', self.prositeProfiles_hits))
                writer.writerow(('Number of different ProSiteProfiles protein domains found:', len(self.unique_prositeProfiles_hits)))
                writer.writerow('')
                writer.writerow(('Total number of ProSitePatterns protein domain families:', self.prositePatterns_hits))
                writer.writerow(('Number of different ProSitePatterns protein domains found:', len(self.unique_prositePatterns_hits)))
                writer.writerow('')
                writer.writerow(('Total number of SUPERFAMILY structural annotations:', self.superfamily_hits))
                writer.writerow(('Number of different SUPERFAMILY structural annotations:', len(self.unique_superfamily_hits)))
                writer.writerow('')
                writer.writerow(('Total number of PRINTS motifs:', self.prints_hits))
                writer.writerow(('Number of different PRINTS motifs found:', len(self.unique_prints_hits)))
                writer.writerow('')
                writer.writerow(('Total number of PANTHER protein families:', self.panther_hits))
                writer.writerow(('Number of different PATNTHER protein families found:', len(self.unique_panther_hits)))
                writer.writerow('')
                writer.writerow(('Total number of Gene3D domains:', self.gene3D_hits))
                writer.writerow(('Number of different Gene3D domains found:', len(self.unique_gene3D_hits)))
                writer.writerow('')
                writer.writerow(('Total number of PIRSF domains:', self.pirsf_hits))
                writer.writerow(('Number of different PIRSF domains found:', len(self.unique_pirsf_hits)))
                writer.writerow('')
                writer.writerow(('Total number of Pfam domains:', self.pfam_hits))
                writer.writerow(('Number of different Pfam domains found:', len(self.unique_pfam_hits)))
                writer.writerow('')
                writer.writerow(('Total number of TIGRFAM domains:', self.tigrfam_hits))
                writer.writerow(('Number of different TIGRFAM domains found:', len(self.unique_tigrfam_hits)))
                writer.writerow('')
                writer.writerow(('Total number of genes with potential Coiled-coil conformation:', self.coils_hits))





parser = argparse.ArgumentParser(description='Script to process and produce some stats about an interproscan run.')
parser.add_argument(dest='interpro_file', metavar='interproFile', nargs=1,help='Interproscan output file in tsv format to be processed.')
parser.add_argument(dest='output_file', metavar='outputFile', nargs=1, help='Output file where some stats will be written.')
args = parser.parse_args()


if __name__ == "__main__":

    hashMap = processInterproRun(args.interpro_file[0])
    interpro_object = InterproDomains()
    interpro_object.analyseHash(hashMap)
    interpro_object.writeReport(args.output_file[0])



