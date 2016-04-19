__author__ = 'pedro'

import argparse
import os
import csv
import numpy as np
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
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
       self.gene_go_terms = defaultdict(set)

       self.unique_kegg_pathways = set()
       self.gene_kegg_data = defaultdict(set)
       self.unique_kegg_enzymes = set()


       self.unique_metacyc_pathways = set()
       self.unique_reactome = set ()
       self.unique_unipathway = set()

       self.average_domains_gene = []

    def analyseHash(self,inputDict):


        self.genes_with_domains = len(inputDict)
        for gene,hits in iter(inputDict.items()):
            self.average_domains_gene.append(len(hits)) #number of hits for this gene


            #each hit has the following fields: analysis[0],signature_id[1], signature description[2], score[3], interpro_id[4], interpro_description[5], go_terms[6]
            #pathway_analysis[7]

            for hit in hits:
                self.total_number_domains +=1


                if len(hit) > 4 and len(hit) < 7: #with interpro annotation
                    self.interpro_hits +=1
                    self.unique_interpro_id.add(hit[4]) #interpro id

                elif len(hit) == 7: #with go annotation
                    if hit[4]:
                        self.interpro_hits +=1
                        self.unique_interpro_id.add(hit[4]) #interpro id


                    if "|" in hit[6]: #gene_go dict
                        gos = hit[6].split("|")
                        for go in gos:
                            self.gene_go_terms[gene].add(go)
                            self.unique_go_terms.add(go) #go id
                    elif hit[6]:
                        self.gene_go_terms[gene].add(hit[6])
                        self.unique_go_terms.add(hit[6]) #go id
                    self.go_terms += 1


                elif len(hit) > 7 : #with pathway annotation
                    if hit[4]:
                        self.interpro_hits +=1
                        self.unique_interpro_id.add(hit[4]) #interpro id


                    if hit[6]:

                        self.go_terms += 1
                        if "|" in hit[6]: #gene_go_dict
                            gos = hit[6].split("|")
                            for go in gos:
                                self.gene_go_terms[gene].add(go)
                                self.unique_go_terms.add(go) #go id
                        else:
                            self.unique_go_terms.add(hit[6])
                            self.gene_go_terms[gene].add(hit[6])


                    self.pathway_annotation +=1
                    ##parse pathway field##
                    #print (hit[7])
                    diff_info = hit[7].split("|")
                    for analysis in diff_info:

                        if "KEGG:" in analysis:
                            #pathway
                            kegg=analysis[5:].strip()
                            self.kegg_pathways += 1
                            kegg_id = kegg.split('+')[0]
                            self.unique_kegg_pathways.add(kegg_id)
                            self.gene_kegg_data[gene].add(kegg)

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


    def writeReport(self,outputBase):
            if os.path.exists(outputBase + "-stats.txt"):
                os.remove(outputBase + "-stats.txt")
            with open(outputBase + "-stats.txt", "w") as csvfile:
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
                writer.writerow(('Number of genes with GO terms found:', len(self.gene_go_terms)))
                writer.writerow('')
                writer.writerow(('Total number of matches with pathway analysis:', self.pathway_annotation))
                writer.writerow(('Number of genes  with KEGG pathway information:', len(self.gene_kegg_data)))
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



    def goMapping(self,outputBase,bio_process,mol_function,cell_component):
        bio_process_hash = {}
        mol_function_hash = {}
        cell_component_hash = {}
        all_gos = []
        with open(bio_process, 'r') as bioFile:

            for line in bioFile:
                line.rstrip()
                go_term = line[:10]
                go_description = line[11:]
                bio_process_hash[go_term] = go_description
                if go_term not in all_gos:
                    all_gos.append(go_term)
                else:
                    logging.info("WARNING: GO %s found more than once." % go_term)
        bioFile.close()

        with open(mol_function, 'r') as molFile:
            for line in molFile:
                line.rstrip()
                go_term = line[:10]
                go_description = line[11:]
                mol_function_hash[go_term] = go_description
                if go_term not in all_gos:
                    all_gos.append(go_term)
                else:
                    logging.info("WARNING: GO %s found more than once." % go_term)
        molFile.close()

        with open(cell_component,'r') as cellFile:
            for line in cellFile:
                line.rstrip()
                go_term = line[:10]
                go_description = line[11:]
                cell_component_hash[go_term] = go_description
                if go_term not in all_gos:
                    all_gos.append(go_term)
                else:
                    logging.warning("WARNING: GO %s found more than once." % go_term)
            cellFile.close()

        with open(outputBase + "-GO-BiologicalProcess.txt", "w") as bp_file:
            with open(outputBase + "-GO-MolecularFunction.txt", "w") as mf_file:
                with open(outputBase + "-GO-CelularComponent.txt", "w") as cc_file:
                    writer_bp = csv.writer(bp_file,dialect=csv.excel_tab)
                    writer_mf = csv.writer(mf_file,dialect=csv.excel_tab)
                    writer_cc = csv.writer(cc_file,dialect=csv.excel_tab)

                    writer_bp.writerow(("#gene_id","#GO_term", "#Description"))
                    writer_mf.writerow(("#gene_id","#GO_term", "#Description"))
                    writer_cc.writerow(("#gene_id","#GO_term", "#Description"))

                    for gene,gos in iter(self.gene_go_terms.items()):
                        for go in gos:
                            if go in bio_process_hash:
                                writer_bp.writerow((gene, go, bio_process_hash[go].rstrip()))
                            elif go in mol_function_hash:
                                writer_mf.writerow((gene, go, mol_function_hash[go].rstrip()))
                            elif go in cell_component_hash:
                                writer_cc.writerow((gene, go, cell_component_hash[go].rstrip()))



    def keggMapping(self,outputBase,keggPathNames):
        kegg_mapping_hash = {}

        with open(keggPathNames, 'r') as keggFile:
            for line in keggFile:
                line.rstrip()
                kegg_ids = line[:5]
                kegg_path_names = line[7:]
                kegg_mapping_hash[kegg_ids] = kegg_path_names.rstrip()
        keggFile.close()

        with open(outputBase + "-KEGG-gene2pathway.txt", "w") as kegg_file:
            writer_kegg = csv.writer(kegg_file,dialect=csv.excel_tab)

            writer_kegg.writerow(("#gene_id","#pathway_name","#pathway_id","#ec_number"))
            for gene, keggData in iter(self.gene_kegg_data.items()):
                for kegg in keggData:
                    data = kegg.split("+")
                    pathway_id = data[0]
                    description = ""
                    if pathway_id in kegg_mapping_hash:
                        description = kegg_mapping_hash[pathway_id]
                    else:
                        logging.warning("Pathway id %s not found in kegg mapping file provided" % pathway_id)


                    for i in range(1,len(data),1):
                        ec_number = data[i]
                        writer_kegg.writerow((gene,description,pathway_id,ec_number))


parser = argparse.ArgumentParser(description='Script to process and produce some stats about an interproscan run. If desirable, some additional analysis can be executed, particularly the KEGG and Gene Ontology mapping.')
parser.add_argument(dest='interpro_file', metavar='interproFile', nargs=1,help='Interproscan output file in tsv format to be processed.')
parser.add_argument(dest='output_basename', metavar='outputBasename', nargs=1, help='Basename for the output files.')
parser.add_argument('-G', '--goMapping', action='store_true', help='Flag to output additional file mapping GO terms found by interpro to categories.')
parser.add_argument('-P', '--pathwayMapping', action='store_true', help='Flag to output additional file mapping genes to KEGG pathways.')
parser.add_argument('-b', '--biologicalProcess', nargs=1, help='Required if -G, file mapping GOs to the Biological Process category.')
parser.add_argument('-m', '--molecularFunction', nargs=1, help='Required if -G, file mapping GOs to the Molecular Function category.')
parser.add_argument('-c', '--cellularCompoment',nargs=1, help='Required if -G, file mapping GOs to the Cellular Component category.')
parser.add_argument('-n', '--nameKeggPathway', nargs=1, help='Required, if -P is set, file mapping KEGG pathways IDs to their names separated by space.')
args = parser.parse_args()


if __name__ == "__main__":

    if args.goMapping:
        if not args.biologicalProcess or not args.molecularFunction or not args.cellularCompoment:
            logging.error("ERROR: If '-G', please set arguments -b, -m and -c!")
            sys.exit(2)
        elif not os.path.exists(args.biologicalProcess[0]):
            logging.error("ERROR: File %s does not exist.", args.biologicalProcess[0])
            exit(2)
        elif not os.path.exists(args.molecularFunction[0]):
            logging.error("ERROR: File %s does not exist.", args.molecularFunction[0])
            exit(2)
        elif not os.path.exists(args.cellularCompoment[0]):
            logging.error("ERROR: File %s does not exist.", args.cellularCompoment[0])
            exit(2)

    elif args.biologicalProcess or args.molecularFunction or args.cellularCompoment:
        logging.warning("ERROR:You provided at least one of the GO auxiliar files, but '-G' flag was not set. Please set '-G' flag.")
        sys.exit(2)

    if args.pathwayMapping:
        if not args.nameKeggPathway:
            logging.error("ERROR: If '-P', please set argument -n!")
            sys.exit(2)
        elif not os.path.exists(args.nameKeggPathway[0]):
            logging.error("ERROR: File %s does not exist.", args.nameKeggPathway[0])
            sys.exit(2)

    elif args.nameKeggPathway:
        logging.error("ERROR:You provided KEGG file mapping pathay IDs to their name, but '-P' flag was not set. Please set '-P' flag.")
        sys.exit(2)


    logging.info("Processing interpro output file..")
    hashMap = processInterproRun(args.interpro_file[0])
    interpro_object = InterproDomains()
    interpro_object.analyseHash(hashMap)

    if args.goMapping:
        logging.info("Processing GO mapping files..")
        interpro_object.goMapping(args.output_basename[0],args.biologicalProcess[0],args.molecularFunction[0],args.cellularCompoment[0])

    if args.pathwayMapping:
          logging.info("Processing KEGG mapping files..")
          interpro_object.keggMapping(args.output_basename[0],args.nameKeggPathway[0])

    logging.info("Writing final report.")
    interpro_object.writeReport(args.output_basename[0])