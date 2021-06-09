import argparse
from _collections import defaultdict
import plotnine as p9
import pandas as pd
import numpy as np
import os
from typing import Union


def process_vastools_proportions(tidy_files: list,
                                 event_ids: list,
                                 comparisons: list,
                                 groups: list = None):
    """
    :param list tidy_files: Tidy files, one for each comparison
    :param list event_ids: List of events ids
    :param list comparisons: Names for each comparison
    :param list groups: Ordered list of groups that were compared
    in the vastools tidy run
    :return list: Data to plot
    """
    out_data = []
    group1, group2 = "", ""

    for idx_file, file in enumerate(tidy_files):

        f = open(file, "r")

        for idx_line, line in enumerate(f):

            fields = line.rstrip().split("\t")
            if idx_line == 0:
                if groups:
                    assert all([x in [fields[1], fields[2]] for x in groups]), "Groups provided not found " \
                                                                               "in tidy file"

                    group1 = fields[fields.index(groups[0])]
                    group2 = fields[fields.index(groups[1])]
                    idx_g1 = fields.index(groups[0])
                    idx_g2 = fields.index(groups[1])
                else:
                    group1 = fields[1]
                    group2 = fields[2]
                    idx_g1 = 1
                    idx_g2 = 2

                continue

            gene_name = fields[0].split("=")[0]

            id = fields[0].split("=")[1]
            if id in event_ids:
                comparison = comparisons[idx_file]

                inc_g1 = float(fields[idx_g1])
                inc_g2 = float(fields[idx_g2])

                out_data.append([id, comparison, group1, inc_g1,
                                 group2, inc_g2, gene_name])

    return out_data


def process_rmats_proportions(rmats_file: str,
                              event_ids: list,
                              groups: list):
    """
    Process rMATs file

    :param str rmats_file: Differential spliced events
    processed by the process_rmats_maser.py script
    :param str event_ids: Event IDs to plot
    :param list groups: Groups compared in a single rmats run
    """
    out_data = []
    f = open(rmats_file, 'r')

    for line in f:

        if line.split("\t")[0] in event_ids:
            fields = line.rstrip().split("\t")
            id = fields[0]
            fdr = float(fields[2])
            gene_id = fields[3]
            gene_name = fields[4]
            ev_type = fields[5]
            comparison = "single_group" if not fields[6] else fields[6]
            inc_g1 = round(np.mean([float(x) for x in fields[7].split(",")]), 3)
            inc_g2 = round(np.mean([float(x) for x in fields[8].split(",")]), 3)

            coords = ':'.join(id.split(":")[2:])
            chrom = id.split(":")[0]
            strand = id.split(":")[1]

            out_data.append([id, comparison, groups[0], inc_g1, fdr,
                             groups[1], inc_g2, gene_id, coords,
                             gene_name, chrom, strand, ev_type])
    return out_data


def get_relevant_type_majiq(lsvtype):
    """
    Get junction types within a LSV
    :param list lsvtype: List with boolean
    columns for specific splicing junctions
    :return list: Types of junctions present
    in a given LSV
    """
    d = {0: "A5SS", 1: "A3SS", 2: "ES"}
    return [d[i] for i, v in enumerate(lsvtype) if v == 'True']


def process_majiq_deltapsi(deltapsi_files: list,
                           event_ids: list,
                           comparisons: list,
                           voila_file: str = None):
    """
    Process LSVs in order to extract splicing binary
    information to be plotted.

    :param list deltapsi_files: List of deltapsi files
    :param list event_ids: List of event ids
    :param list comparisons: Names for each comparison
    :param str voila_file: Voila file to extract
    exon coordinates (and gene names). Not mandatory
    """

    out_data = defaultdict(list)
    for idx_file, file in enumerate(deltapsi_files):

        f = open(file, "r")
        group1, group2 = "", ""
        for idx_line, line in enumerate(f):
            fields = line.rstrip().split("\t")

            lsv_id = fields[1]

            if idx_line == 0:
                group1 = fields[6].split()[0]
                group2 = fields[7].split()[0]

            if lsv_id in event_ids:
                lsv_type_str = ';'.join(get_relevant_type_majiq(fields[8:11]))
                lsv_source_or_target = fields[2]
                gene_id = fields[0].split('.')[0]
                if len(fields) == 15:
                    lsv_type_str += "RI" if lsv_type_str == "" else ";RI"

                # starts by (s or t) being source or target
                # each '|' is a new junction representation and if there is intron_retention the last character is 'i'
                # each junction is represented by  XeY.ZoK where
                # X is the ordinal splice site in the reference exon
                # Y is the ordinal exon connecting the lsv
                # Z is the ordinal splice site in exon Y
                # K is the total number of splice sites that Y has

                # If multiple files, we will just keep the same junction between comparisons
                if lsv_id in out_data.keys():
                    dpsi_relevant_idx = out_data[lsv_id][0][-1]

                else:
                    # dPSI magnitudes refer to 2nd group vs 1st group
                    dpsi_relevant_idx = [i for i, x in enumerate(fields[3].split(";")) if abs(float(x)) > 0.2]

                psi_group1 = [float(x) for i, x in enumerate(fields[6].split(";"))]  # if i in dpsi_relevant_idx]
                psi_group2 = [float(x) for i, x in enumerate(fields[7].split(";"))]  # if i in dpsi_relevant_idx]
                prob_g1 = [float(x) for i, x in enumerate(fields[4].split(";"))]  # if i in dpsi_relevant_idx]
                prob_g2 = [float(x) for i, x in enumerate(fields[4].split(";"))]  # if i in dpsi_relevant_idx]
                jx_coord = [x for i, x in enumerate(fields[13].split(";"))]  # if i in dpsi_relevant_idx]

                assert len(psi_group1) == len(psi_group2), "Wrong number of jx " \
                                                           "between groups for {} jx".format(lsv_id)

                # target exon - most left junction refers to the skipping
                if lsv_source_or_target.split("|")[0] == "t":

                    if len(psi_group1) == 2:
                        inc_g1 = round(psi_group1[1] * 100, 2)
                        inc_g2 = round(psi_group2[1] * 100, 2)
                        prob_inc_g1 = prob_g1[1]
                        prob_inc_g2 = prob_g2[1]
                        # exc_g1 = psi_group1[0]
                        # exc_g2 = psi_group2[0]

                    else:
                        inc_g1 = [round(x * 100, 2) for x in psi_group1[1:]]
                        inc_g2 = [round(x * 100, 2) for x in psi_group2[1:]]
                        prob_inc_g1 = [x for x in prob_g1[1:]]
                        prob_inc_g2 = [x for x in prob_g2[1:]]

                        # exc_g1 = psi_group2[0]
                        # exc_g2 = psi_group2[0]

                # source exon - most left junction refers to the inclusion
                elif lsv_source_or_target.split("|")[0] == "s":

                    if len(psi_group1) == 2:
                        inc_g1 = round(psi_group1[0] * 100, 2)
                        inc_g2 = round(psi_group2[0] * 100, 2)
                        prob_inc_g1 = prob_g1[0]
                        prob_inc_g2 = prob_g2[0]
                        # exc_g1 = psi_group1[1]
                        # exc_g2 = psi_group2[1]

                    else:
                        inc_g1 = [round(x * 100, 2) for x in psi_group1[:-1]]
                        inc_g2 = [round(x * 100, 2) for x in psi_group2[:-1]]
                        prob_inc_g1 = [x for x in prob_g1[:-1]]
                        prob_inc_g2 = [x for x in prob_g2[:-1]]
                        # exc_g1 = psi_group1[-1]
                        # exc_g2 = psi_group2[-1]

                out_data[lsv_id].append([comparisons[idx_file], group1, inc_g1, prob_inc_g1,
                                         group2, inc_g2, prob_inc_g2, gene_id, ';'.join(jx_coord), dpsi_relevant_idx])

    # Add gene name, chromosome, exon coords
    if voila_file:
        f = open(voila_file, 'r')
        for i, line in enumerate(f):
            if not line.startswith("#"):

                fields = line.rstrip().split("\t")
                if fields[2] in out_data.keys():

                    for i, comparison in enumerate(out_data[fields[2]]):
                        out_data[fields[2]][i] = out_data[fields[2]][i][:-1]
                        out_data[fields[2]][i].extend([fields[0], fields[12], fields[13], fields[15]])

    else:
        out_data = {k: v[:-1] for k, v in out_data.items()}

    return out_data


def _draw(df: pd.DataFrame,
          figsize: tuple,
          textsize: int,
          textangle: int) -> p9.ggplot:
    """
    Draw plot
    :param pd.DataFrame df: Data to plot
    :param dict colors: Colors to represent inclusion and exclusion values
    """

    colors = {'inc': 'skyblue',
              'inc_1': 'skyblue', 'inc_2': 'brown',
              'inc_3': 'darkgrey', 'inc_4': 'darkblue',
              'exc': 'mistyrose'}

    p1 = (p9.ggplot(df) +
          p9.geom_col(color='black', size=0.5, mapping=p9.aes(x='class', y='psi', fill='group')) +
          p9.facet_grid('comparison ~ gene_name') +
          p9.scale_fill_manual(values=colors) +
          p9.theme(figure_size=figsize,
                   text=p9.element_text(size=textsize),
                   legend_position="none",
                   axis_text_x=p9.element_text(angle=textangle)) +

          p9.xlab('') +
          p9.ylab('PSI') +
          p9.ylim(0, 100))

    return p1


def plot_psi_proportions(data: Union[dict, list], outdir: str, tool: str):
    """
    Draw
    :param Union[dict, list data: Dict with all the info needed to plot
    :param str outdir: Output directory
    :param str tool: Tool for which the data refers to
    """
    def _update_majiq_inc_group(df: pd.DataFrame):
        """
        Updates majiq inclusion group when there are multiple jx supporting it
        :param pd.DataFrame df: PSI values for a single group
        :return pd.DataFrame: Updated df
        """

        if df.shape[0] > 1 and df.iloc[0].group == "inc":
            df['group'] = df.group.str.cat(df.groupby('group').cumcount().add(1).astype(str), sep="_")

        return df

    def _format_df(long_df: pd.DataFrame) -> pd.DataFrame:
        """
        Formats some columns for proper plotting
        :param pd.DataFrame long_df: Melted df with PSI values
        :return pd.DataFrame: Fixed df
        """
        long_df[['group', 'class']] = long_df.ev_type.str.split('_', expand=True)
        # SPECIFIC FOR CBE DATA
        long_df['g2'] = long_df['g2'].str.replace('_B6', '').str.replace('_BC', '')

        group_names = {'g1': long_df.iloc[0]['g1'],
                       'g2': long_df.iloc[0]['g2']}

        long_df['class'] = np.where(long_df['class'] == 'g1', group_names['g1'], group_names['g2'])
        long_df['class'] = pd.Categorical(long_df['class'],
                                          categories=list(group_names.values()),
                                          ordered=True)
        return long_df

    if tool == "vastools":
        df = pd.DataFrame.from_records(data, columns=['id', 'comparison', 'g1',
                                                      'inc_g1', 'g2', 'inc_g2',
                                                      'gene_name'])

        df['exc_g1'] = 100 - df['inc_g1']
        df['exc_g2'] = 100 - df['inc_g2']

        long_df = pd.melt(df, id_vars=['id', 'gene_name', 'comparison', 'g1', 'g2'],
                          value_vars=['inc_g1', 'inc_g2', 'exc_g1', 'exc_g2'],
                          var_name='ev_type',
                          value_name='psi')

        long_df = _format_df(long_df)

    if tool == "rmats":
        df = pd.DataFrame.from_records(data, columns=['id', 'comparison', 'g1',
                                                      'inc_g1', 'fdr', 'g2',
                                                      'inc_g2', 'gene_id', 'exon_coord',
                                                      'gene_name', 'chr', 'strand', 'event_type'])

        df[['inc_g1', 'inc_g2']] = df[['inc_g1', 'inc_g2']].apply(lambda x: x * 100)

        df['exc_g1'] = 100 - df['inc_g1']
        df['exc_g2'] = 100 - df['inc_g2']

        long_df = pd.melt(df, id_vars=['id', 'gene_id', 'gene_name', 'comparison', 'g1', 'g2'],
                          value_vars=['inc_g1', 'inc_g2', 'exc_g1', 'exc_g2'],
                          var_name='ev_type',
                          value_name='psi')

        long_df = _format_df(long_df)

    if tool == "majiq":
        to_df = []

        for id, majiq_info in data.items():
            for comparison in majiq_info:
                to_df.append([id] + comparison)

        df = pd.DataFrame.from_records(to_df)

        if len(list(df)) == 14:
            df.columns = ['id', 'comparison', 'g1', 'inc_g1',
                          'prob_g1', 'g2', 'inc_g2', 'prob_g2',
                          'gene_id', 'jx_coords', 'gene_name',
                          'chr', 'strand', 'exon_coord']
        else:
            df.columns = ['id', 'comparison', 'g1', 'inc_g1',
                          'prob_g1', 'g2', 'inc_g2', 'prob_g2',
                          'gene_id', 'jx_coords']

        df['exc_g1'] = df.inc_g1.apply(lambda x: 100 - x if isinstance(x, float) else round(100 - np.sum(x), 2))
        df['exc_g2'] = df.inc_g2.apply(lambda x: 100 - x if isinstance(x, float) else round(100 - np.sum(x), 2))

        to_explode_cols = ['inc_g1', 'inc_g2']
        df = df.drop(columns=['prob_g1', 'prob_g2'])

        df = df.set_index([col for col in list(df) if col not in to_explode_cols]).apply(
            pd.Series.explode).reset_index()

        long_df = pd.melt(df, id_vars=['id', 'gene_id', 'gene_name', 'comparison', 'g1', 'g2'],
                          value_vars=['inc_g1', 'inc_g2', 'exc_g1', 'exc_g2'],
                          var_name='ev_type',
                          value_name='psi')

        long_df = _format_df(long_df)

        long_df = long_df.groupby(['id', 'comparison', 'ev_type']).apply(lambda x: _update_majiq_inc_group(x))

    long_df.sort_values(by=['gene_name']).to_csv(os.path.join(outdir, 'data.tsv'),
                                               sep="\t",
                                               index=False)

    long_df['psi'] = pd.to_numeric(long_df.psi)
    p1 = _draw(long_df, figsize=(11, 4), textsize=10, textangle=60)
    p1.save(os.path.join(outdir, 'all.pdf'), verbose=False)

    for name, group in long_df.groupby('id'):
        p1 = _draw(group, figsize=(2, 2), textsize=10, textangle=0)
        p1.save(os.path.join(outdir, name + ".pdf"), verbose=False)


def main():
    parser = argparse.ArgumentParser(description='Script to plot splicing ratios of individual '
                                                 'events provided by event IDs.')

    parser.add_argument(dest='event_ids', help='Path to a file with event IDs to be individually plotted')
    parser.add_argument('-o', '--output_dir', required=True, help="Output dir")
    # parser.add_argument('-individual', action='store_true', help='Whether plots should be drawn separately '
    #                                                              'for each event ID. Default, stack all plots '
    #                                                              'tohether on a single image')

    parser.add_argument("-r", "--rmats", help='Path to the output file of process_rmats_maser script '
                                              'called *psi_proportions.tsv, where all the inclusion '
                                              'levels for all comparisons are displayed, regardless of '
                                              'the significance.')

    parser.add_argument("-g", "--groups_rmats", metavar="", nargs=2,
                        help='Groups that were compared in the rMATs run. If multiple comparisons '
                             'comparisons were done, the same groups will apply.')

    parser.add_argument('-v', '--vastools', nargs="+", metavar="",
                        help='Path to the output file of vastools tidy '
                             'where inclusion levels for each group are '
                             'written.')

    parser.add_argument('--names_vastools', nargs='+', metavar='', help='Names for the comparisons represented on '
                                                                        'each tidy file provided '
                                                                        'by the vastools tool.')

    parser.add_argument('--groups_vastools', nargs=2, metavar='', help='Ordered groups that were compared in '
                                                                       'the vastools tidy run. Although '
                                                                       'tidy has that info, it may not be '
                                                                       'ordered as desired.')

    parser.add_argument("-m", "--majiq", nargs="+", help='Path to the output file of majiq deltapsi. Multiple '
                                                         'files can be provided to represent multiple deltapsi '
                                                         'comparisons.')

    parser.add_argument('--majiq_voila', help="Path to the output of a single majiq voila. This is only "
                                              "used to extract exon coordinates, it's not mandatory. If "
                                              "multiple comparisons are going to be plotted, it assumes "
                                              "that a single majiq build step including all the samples "
                                              "was done, therefore ensuring the LSV IDs represent the "
                                              "same event across multiple comparisonOrder")

    parser.add_argument('--names_majiq', nargs='+', metavar='', help='Names for the comparisons represented on '
                                                                     'each differential splicing file provided '
                                                                     'by the majiq tool.')

    args = parser.parse_args()

    if all(x is None for x in [args.vastools, args.majiq, args.rmats]):
        raise ValueError('Please provide input files for the tool output the IDs provided.')

    if len([x for x in [args.vastools, args.majiq, args.rmats] if x]) != 1:
        raise ValueError('Only files from tool must be provided for the given event IDs.')

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    id_list = [line.rstrip() for line in open(args.event_ids)]

    if args.majiq:
        if not args.names_majiq:
            raise ValueError('Names for the comparisons (--names_majiq) must be given when --majiq is set.')

        assert len(args.majiq) == len(args.names_majiq), "Names of the comparisons must be " \
                                                         "given to each majiq file provided."

        data_to_plot = process_majiq_deltapsi(args.majiq, id_list, args.names_majiq, args.majiq_voila)
        plot_psi_proportions(data_to_plot, args.output_dir, 'majiq')

    if args.rmats:
        if not args.groups_rmats:
            raise ValueError('Names for the groups (--groups_rmats) must be given when --rmats is set.')

        data_to_plot = process_rmats_proportions(args.rmats, id_list, args.groups_rmats)

        plot_psi_proportions(data_to_plot, args.output_dir, 'rmats')

    if args.vastools:
        if not args.names_vastools:
            raise ValueError('Names for the groups (--names_vastools) must be given when --vastools is set.')

        assert len(args.vastools) == len(args.names_vastools), "Names of the comparisons must be " \
                                                               "given to each vastools tidy file provided."

        if args.groups_vastools:
            assert len(args.groups_vastools) == 2, "Two groups are required when using --group_vastools arg"

        data_to_plot = process_vastools_proportions(args.vastools, id_list, args.names_vastools, args.groups_vastools)
        plot_psi_proportions(data_to_plot, args.output_dir, 'vastools')


if __name__ == "__main__":
    main()
