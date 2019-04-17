import argparse
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from gtfparse import read_gtf
from pyensembl import Genome
import binascii
import gzip
import numpy as np
import pandas as pd
import re
from collections import defaultdict

def extend_attributes(x, cols):
    y=x.replace("\"", "").rstrip()
    R_SEMICOLON = re.compile(r'\s*;\s*')
    fields=[i for i in re.split(R_SEMICOLON, y) if i.strip()]
    attr=defaultdict(list)
    for f in fields:
        attr[re.split('\s+', f)[0]].append(re.split('\s+', f)[1])

    return [';'.join(attr[col]) if col in attr.keys() else np.nan for col in cols]


def get_extend_attributes(df):
    additional_cols=[
        'gene_id',
        'gene_name',
        'gene_type',
        'level',
        'ont',
        'tag',
        'transcript_id',
        'transcript_name',
        'transcript_support_level',
        'transcript_type',
        'exon_id',
        'exon_number'
    ]

    df_new = df['attribute'].apply(extend_attributes, args=(additional_cols,)).to_frame()
    df_with_attributes = pd.DataFrame(df_new['attribute'].values.tolist(), index=df.index, columns=additional_cols)
    df_final=pd.concat([df.drop(['attribute'], axis=1), df_with_attributes], axis=1, sort=False)
    df_final.to_csv("gtf_pandas.csv.gz", compression='gzip')
    return df_final

def parse_gtf(gtf):

    gtf_cols = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    df = pd.read_csv(
        gtf,
        sep="\t",
        comment="#",
        names=gtf_cols,
        skipinitialspace=True,
        skip_blank_lines=True,
        error_bad_lines=True,
        warn_bad_lines=True,
        engine="c",
        dtype={
            "start": np.int64,
            "end": np.int64,
            "score": np.float32,
            "seqname": str,
        },
        na_values=".")

    return get_extend_attributes(df)


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def extract_transcripts(df, biotype):
    """
Check length of removed transcripts
    :param df:
    :param biotype:
    :return:
    """

    if biotype == "protein_coding":
        canonical_tags = ['basic', 'CCDS', 'MANE_Select', 'appris_principal_1', 'appris_principal_2', 'appris_principal_3',
                          'appris_principal_4', 'appris_principal_5']
    else:
        canonical_tags = ['basic']

    genes_biotype_filtered = df.loc[(df['gene_type']) == biotype]
    grouped_genes = genes_biotype_filtered.groupby('gene_id')
    transcripts_per_gene = grouped_genes.apply(lambda x: x[x['feature'] == 'transcript'])

    canonical = pd.DataFrame(columns=list(df))
    for gene in transcripts_per_gene.index.get_level_values('gene_id').unique():
        logging.info("Gene {}".format(gene))
        df_g = transcripts_per_gene.ix[gene]
        candidate=[]
        for index, transcript in df_g.loc[(df['transcript_type']) == biotype].iterrows():

            if all(x in canonical_tags for x in transcript['tag'].split(";")):
                candidate.append(transcript)

        if len(candidate) == 0:
            logging.info("Problem here")
            print(df_g['tag'])

        elif len(candidate) > 1:
            #index with the highest number of tags
            idx = max([(i, len(cand['tag'].split(";"))) for i, cand in enumerate(candidate)], key=lambda item: item[1])[0]
            canonical = canonical.append(candidate[idx], ignore_index=True)

        elif len(candidate) == 1:
            canonical = canonical.append(candidate, ignore_index=True)

def main():
    parser = argparse.ArgumentParser(description='Utility to extract canonical transcript information from an input GTF file')
    parser.add_argument('--gtf_file', metavar='gtf', help='Annotation file to be analyzed.')
    parser.add_argument('--pandas_db', metavar='pandas_db', help='Pandas dataframe of the gtf file. Useful to skip initial processing of attributes field')
    parser.add_argument('--biotype', metavar='biotype', default="protein_coding", help='Extract transcripts belonging to genes of the given biotype.')

    args = parser.parse_args()
    if args.pandas_db:
        logging.info("Reading pandas db file")
        df = pd.read_csv(args.pandas_db)
    else:
        logging.info("Creating pandas db from gtf")
        df = parse_gtf(args.gtf_file)

    extract_transcripts(df, args.biotype)

if __name__ == "__main__":
    main()