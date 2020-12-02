# Script that processes the output of bcftools split-vep plugin
# in which multiple consequences for the same variant are split
# with ','. Example of a valid command to produce a valid file
# for this script:

# bcftools +split-vep -a ANN -f "%CHROM\t%POS\t%REF\t%ALT\t%Consequence
# \t%SYMBOL\t%Gene\t%HGVSc\t%gnomADg_AF\t%existing_InFrame_oORFs\t%existing_OutOfFrame_oORFs
# \t%existing_uORFs\t%five_prime_UTR_variant_annotation\t%five_prime_UTR_variant_consequence\n"
# 2_with_UTRannotations.vcf.gz > 3_UTRannot_split.tsv

import argparse
import pandas as pd
import numpy as np


def process_genotypes(df: pd.DataFrame, genotypes: str):
    """
    Process genotypes and merge sample info into df

    :param pd.DataFrame df: UTRannotator processed df
    :param genotypes: File with a map between the variants
        in df and corresponding genotypes
    :return pd.DataFrame: Final df to write output file
    """
    df_genotypes = pd.read_csv(genotypes, sep="\t")
    return pd.merge(df, df_genotypes, how='left', on=['POS', 'REF', 'ALT'])


def process_file(filename: str):
    """
    Processes input file and returns only non-null
    consequences from the UTRannotator VEP plugin

    :param str filename: Input file
    :return pd.DataFrame
    """
    _file = open(filename, 'r')
    static_field_names = ['#CHROM', 'POS', 'REF', 'ALT', 'AF_in_cohort', 'gnomADg_AF']
    dynamic_fields_names = []
    fields_map_idx = {}
    target_fields = ['existing_InFrame_oORFs', 'existing_OutOfFrame_oORFs',
                     'existing_uORFs', 'five_prime_UTR_variant_annotation',
                     'five_prime_UTR_variant_consequence']

    final_df = pd.DataFrame()
    for line in _file:
        line = line.rstrip()
        if line.startswith('#'):
            fields = line.split()
            fields_map_idx = {f: i for i, f in enumerate(fields)}
            idx_map_fields = {y: x for x, y in fields_map_idx.items()}
            static_idx = [v for k, v in fields_map_idx.items() if k in static_field_names]
            conseq_idx = fields_map_idx['Consequence']

        else:
            fields = line.split()
            static_fields_values = [fields[i] for i in static_idx]

            utr5_consequences_idx = [i for i, j in enumerate(fields[conseq_idx].split(','))
                                    if '5_prime_UTR_variant' in j]
            
            variant_out = []
            for i, _f in enumerate(fields):
                # if it's a field from VEP
                if i not in static_idx:

                    if idx_map_fields[i] not in dynamic_fields_names:
                        dynamic_fields_names.append(idx_map_fields[i])
                    _value = [v for _i, v in enumerate(_f.split(',')) if _i in utr5_consequences_idx]
                    variant_out.append(_value)

            n_utr5_conseq = len(variant_out[0])
            _static_exploded = pd.DataFrame([static_fields_values] * n_utr5_conseq)
            _to_concat = pd.DataFrame(variant_out).T
            final_df = final_df.append(other=pd.concat([_static_exploded, _to_concat], axis=1))

    final_df.columns = static_field_names + dynamic_fields_names
    return final_df


def main():
    parser = argparse.ArgumentParser(description='Script to split per-transcript annotations of 5\'UTR '
                                                 'variants produced by UTRannotator.')
    parser.add_argument(dest='infile', help='Path to the tab delimited file produced by '
                                            'bcftools split-vep plugin')
    parser.add_argument('-g', '--genotypes', help='File with genotypes information for the variants passed')
    args = parser.parse_args()
    df = process_file(args.infile)

    if args.genotypes:
        df["AF_in_cohort"] = pd.to_numeric(df["AF_in_cohort"])
        df["POS"] = df['POS'].astype(np.int64)
        df = df[df.AF_in_cohort < 0.01]
        df = process_genotypes(df, args.genotypes)

    df.to_csv('test.tsv', sep="\t", index=False)


if __name__ == "__main__":
    main()

