import argparse
import os
import pandas as pd


def process_files(list_files: list):
    final_dict = {}
    for file in list_files:

        name = os.path.basename(os.path.splitext(file)[0])
        d = pd.read_excel(file).set_index('gene_id').to_dict(orient='index')
        for gene_id, values in d.items():
            log2fc = values['log2FoldChange']
            comparison_info = '{} (log2FC {})'.format(name, str(round(log2fc, 3)))
            if gene_id in final_dict.keys():
                _tmp_val = final_dict[gene_id]
                _tmp_val.append(comparison_info)
                final_dict[gene_id] = _tmp_val
            else:
                symbol = values['symbol']
                description = values['description']
                final_dict[gene_id] = [symbol, description, comparison_info]
    return final_dict


def main():
    parser = argparse.ArgumentParser(description='Script to check the interception of the'
                                                 ' differential expressed features present '
                                                 'in different pairwise tests.')
    parser.add_argument(dest='input_files', metavar='diffExp_files', nargs='+',
                        help='Excel files representing a given comparison. 1st column will be used '
                             'as the feature ID (e.g. ensembl geneID). Expected format of '
                             'the output files: gene_id, symbol, description, log2FC, pvalue, padj')
    parser.add_argument('-o', '--output', required=True, help="Output tsv file")
    args = parser.parse_args()

    if len(args.input_files) < 2:
        raise ValueError('You should specify at least two files to compare.')

    final = process_files(args.input_files)
    out = open(args.output, 'w')
    out.write('#gene_id\tsymbol\tdescription\tcomparisons_with_gene\n')
    for gene_id, values in final.items():
        out.write(gene_id + "\t" + '\t'.join(values) + '\n')
    out.close()


if __name__ == "__main__":
    main()
