import argparse
import os
from cyvcf2 import VCF, Writer
import pysam
import mygene
import numpy as np
import sys, re


def checkIndex(expectofile):
    if not os.path.isfile(expectofile + '.tbi'):
        print('Index file is required. Make sure that {} file exists in the same path'.format(expectofile + '.tbi'))
        exit(1)


def add_fields_to_header(vcf_data, directionality):

    vcf_data.add_info_to_header({'ID': 'ExPecto',
                                 'Description': 'A deep neural network aimed to predict the efffect of variants on gene '
                                                'expression',
                                 'Type': 'Character',
                                 'Number': '.'})

    vcf_data.add_info_to_header({'ID': 'ExPecto_tissue', 'Description': 'Top N tissues for which the variant effect is'
                                                                        ' higher.',
                                 'Type': 'Character',
                                 'Number': '.'})

    vcf_data.add_info_to_header(
        {'ID': 'ExPecto_gene_name',
         'Description': 'Gene name for which ExPecto gave predictions',
         'Type': 'Character',
         'Number': '.'})

    if directionality:

        vcf_data.add_info_to_header(
            {'ID': 'ExPecto_sum',
             'Description': 'Sum of absolute mutation effects within 1KB of a TSS for a gene. ',
             'Type': 'Character',
             'Number': '.'})

        vcf_data.add_info_to_header(
            {'ID': 'ExPecto_pos_constraint',
             'Description': 'Probabilities of positive constraint for a given gene in a given tissue. (high prob ->'
                            ' active expression -> contraint towards higher expression levels)',
             'Type': 'Character',
             'Number': '.'})

        vcf_data.add_info_to_header(
            {'ID': 'ExPecto_neg_constraint',
             'Description': 'Probabilities of positive constraint for a given gene in a given tissue. (high prob ->'
                            ' repressed expression -> contraint towards lower expression levels)',
             'Type': 'Character',
             'Number': '.'})


def file_to_dict(infile, return_header=False):
    d = {}
    nline = 0
    tissue_names = ""
    for line in infile:
        fields = line.rstrip().split("\t")
        if nline == 0 and return_header:
            tissue_names = fields[1:]
        else:
            d[fields[0]] = fields[1:]
        nline += 1

    return d if not return_header else d, tissue_names


def get_variation_direction(gene_name, top_N_tissues, tissues_direct, direct, neg, pos):
    idx_new_files, gene_direct_score, gene_neg_constraint, gene_pos_constraint = [], [], [], []
    for tissue in top_N_tissues:
        if tissue in tissues_direct:
            idx_new_files.append(tissues_direct.index(tissue))
        else:
            idx_new_files.append("-")

    gene_direct_score = [direct[gene_name][i] if i != "-" else "." for i in idx_new_files]
    gene_neg_constraint = [neg[gene_name][i] if i != "-" else "." for i in idx_new_files]
    gene_pos_constraint = [pos[gene_name][i] if i != "-" else "." for i in idx_new_files]

    return gene_direct_score, gene_neg_constraint, gene_pos_constraint, idx_new_files


def processVCF(invcf, out, expecto, topn, directionality, skipVariants):
    expecto_directionality_dir = os.path.join(os.path.dirname(expecto), "variation_potential_scores")
    vcf_data = VCF(invcf, gts012=True)
    tbx_expecto = pysam.TabixFile(expecto)

    add_fields_to_header(vcf_data, directionality)
    w = Writer(out, vcf_data)
    mg = mygene.MyGeneInfo()
    tissue_names_with_spaces = tbx_expecto.header[0].split("\t")[5:]
    tissue_names = [t.replace(" ", "_").rstrip() for t in tissue_names_with_spaces]
    if directionality:
        direct_scores = open("{}/directionality_scores.txt".format(expecto_directionality_dir), "r")
        negative_constraint = open("{}/negative_constraint.probability.txt".format(expecto_directionality_dir), "r")
        positive_constraint = open("{}/positive_constraint.probability.txt".format(expecto_directionality_dir), "r")

        direct, tissues_direct_with_spaces = file_to_dict(direct_scores, return_header=True)
        pos = file_to_dict(negative_constraint)[0]
        neg = file_to_dict(positive_constraint)[0]
        tissues_direct = [t.replace(" ", "_").rstrip() for t in tissues_direct_with_spaces]

        direct_scores.close()
        negative_constraint.close()
        positive_constraint.close()

    for record in vcf_data:
        score_retrieved = False
        try:
            for row in tbx_expecto.fetch(record.CHROM, record.start, record.end):
                fields = row.split()
                position = fields[1]
                ref = fields[2]
                alt = fields[3]
                gene_id = fields[4]
                scores = list(map(float, fields[5:]))

                if int(int(position) == record.POS and ref == record.REF and alt == record.ALT[0]):
                    gene_name = mg.query(gene_id, fields="symbol", as_dataframe=True, size=1, species="human").symbol[0]
                    top_N_val_idx = np.argsort([abs(i) for i in scores])[-topn:]
                    top_N_val = [str(scores[i]) for i in top_N_val_idx][::-1]
                    top_N_tissues = [tissue_names[i] for i in top_N_val_idx][::-1]
                    record.INFO["ExPecto"] = ','.join(top_N_val)
                    record.INFO["ExPecto_tissue"] = ','.join(top_N_tissues)
                    record.INFO["ExPecto_gene_name"] = gene_name
                    if directionality:

                        gene_directionality, neg_constraint_prob, pos_constraint_prob, tissues_new = get_variation_direction(
                            gene_name, top_N_tissues, tissues_direct, direct, neg, pos)

                        record.INFO["ExPecto_sum"] = ','.join(gene_directionality)
                        record.INFO["ExPecto_pos_constraint"] = ','.join(neg_constraint_prob)
                        record.INFO["ExPecto_neg_constraint"] = ','.join(pos_constraint_prob)
                        score_retrieved = True
                    break

        except ValueError:
            print("Something wrong in variant {}, {}, {},{}".format(record.CHROM, record.POS, record.REF, record.ALT))
            # record.INFO["ExPecto"] = "."
            # record.INFO["ExPecto_tissue"] = '.'
            # record.INFO["ExPecto_gene_name"] = '.'
            # if directionality:
            #     record.INFO["ExPecto_sum"] = '.'
            #     record.INFO["ExPecto_pos_constraint"] = '.'
            #     record.INFO["ExPecto_neg_constraint"] = '.'

        if skipVariants and score_retrieved:
            w.write_record(record)
        elif not skipVariants:
            w.write_record(record)


def main():
    parser = argparse.ArgumentParser(description='Script to aid on the task of add ExPecto scores to VCF files.')
    parser.add_argument(dest='vcf', default="-", help='Path to the vcf (uncompressed)')
    parser.add_argument(dest='out', help='Path to the VCF output file')
    parser.add_argument('--path', '--expectoPath', default="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/custom_data/expecto/hg19",
                        help="Path where ExPecto files are located")
    parser.add_argument("-n", "--number", default=3, type=int,
                        help="Top n number of tissues to report highest expression effects."
                             "Default:3")
    parser.add_argument("-d", "--directionality", action="store_true",
                        help="If set, directionality scores will be included")
    parser.add_argument("-s", "--skip", action="store_true", help="Only ouput variants with retrieved ExPecto scores")
    args = parser.parse_args()

    expecto_filename = "expecto_1kb_TSS_v0.3.txt.gz"
    expecto_fullpath = os.path.join(args.path, expecto_filename)
    if os.path.isdir(args.path) and os.path.isfile(expecto_fullpath):
        checkIndex(expecto_fullpath)
        processVCF(args.vcf, args.out, expecto_fullpath, args.number, args.directionality, args.skip)

    else:
        print(
            "Does the ExPecto directory provided exists? Also, make sure the expecto scores (expecto_1kb_TSS_v0.3.txt.gz)"
            " are stored there.")
        exit(1)


if __name__ == "__main__":
    main()
