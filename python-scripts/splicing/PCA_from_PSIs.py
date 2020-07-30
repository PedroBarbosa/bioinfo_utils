import argparse
import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import plotnine as p9


def process_vasttools(file, groups, allowed_nas, paired):
    df = pd.read_csv(file, sep="\t", index_col=0)
    df.dropna(axis=0, thresh=len(df.columns) - allowed_nas, inplace=True)
    plot_PCA(df, groups, paired, "vast-tools")


def process_rMATS(files, groups, allowed_nas, paired):
    inc_level_col = "IncLevel1"
    events = ["SE", "MXE", "A3SS", "A5SS", "RI"]
    event_type_map = {"SE": "ES"}
    results_dict, groups_dict = {}, {}

    samples_names = list(groups.keys())
    with open(files) as f:
        file_list = f.readlines()
    f.close()

    for ev in events:
        print("Processing {} events ".format(event_type_map.get(ev, ev)))
        ev_file = [os.path.basename(e.rstrip()) for e in file_list if ev in os.path.basename(e.rstrip())][0]
        ev_dict = {}
        assert os.path.isfile(ev_file), "{} is not a valid file".format(ev_file)
        f = open(ev_file)
        header = f.readline().split("\t")
        assert inc_level_col in header, "IncLevel1 column (PSI quantification) should be present in " \
                                        "rMATS {} file.".format(ev_file)
        idx_target_column = header.index(inc_level_col)
        next(f)
        for line in f:
            l = line.rstrip().split("\t")
            psi = l[idx_target_column]
            if any(v for v in psi.split(",") if v != "NA"):
                ev_dict["{}_{}".format(event_type_map.get(ev, ev), l[0])] = l[idx_target_column]
        f.close()
        results_dict[ev] = ev_dict

    list_dfs = []
    for ev_type, psi_dict in results_dict.items():
        df = pd.DataFrame.from_dict(psi_dict, orient="index", columns=["psis"])
        df[samples_names] = df['psis'].str.split(",", expand=True).replace("NA", np.nan)
        del df['psis']
        df.dropna(axis=0, thresh=len(df.columns) - allowed_nas, inplace=True)
        print("Number of {} events that will be used in the PCA: {}".format(ev_type, df.shape[0]))
        list_dfs.append(df)
        plot_PCA(df, groups_dict, paired, ev_type, "rMATS")
    concat_df = pd.concat(list_dfs)
    plot_PCA(concat_df, groups_dict, paired, "all", "rMATS")


def plot_PCA(df, groups, paired, name, tool=None):
    """
    Function to plot PCA based on processed df and a groups dict
    """
    print("Number of {} events that will be used in the PCA: {}".format(name, df.shape[0]))
    pca = PCA(n_components=3)
    PCs = pca.fit_transform(df.T)
    print("Amount of variance that the first PCs contain for {} data: {}".format(name, pca.explained_variance_ratio_))

    PC_df = pd.DataFrame(data=PCs,
                         columns=['PC1', 'PC2', 'PC3'],
                         index=df.T.index)

    PC_df = PC_df.merge(pd.DataFrame.from_dict(groups, orient='index', columns=['Group']),
                        left_index=True,
                        right_index=True).rename_axis('Sample').reset_index()

    for pc_pair in [("PC1", "PC2"), ("PC2", "PC3")]:
        if paired:
            PC_df['Ind'] = PC_df['Sample'].apply(lambda x: x.split("_")[0])
            p1 = (p9.ggplot(PC_df, p9.aes(x=pc_pair[0], y=pc_pair[1], color="Ind", shape='Group'))
                  + p9.geom_point(size=6, alpha=0.7, position=p9.position_dodge(width=0.3, preserve="total")))
        else:
            p1 = (p9.ggplot(PC_df, p9.aes(x=pc_pair[0], y=pc_pair[1], fill="Group", shape='Sample'))
                  + p9.geom_point(size=6, alpha=0.7, position=p9.position_dodge(width=0.3, preserve="total")))

        if tool is not None:
            output = "{}_{}_{}vs{}.pdf".format(tool, name, pc_pair[0], pc_pair[1])
        else:
            output = "{}_{}vs{}.pdf".format(name, pc_pair[0], pc_pair[1])
        p1.save(output, verbose=False)


def main():
    parser = argparse.ArgumentParser(description='Script to generate PCAs from inclusion levels of splicing events '
                                                 'quantified by different tools.')
    parser.add_argument(dest='tool', choices=("rmats", "majiq", "vastools"), help="Generate samples PCAs from "
                                                                                  "quantifications of given tool")
    parser.add_argument('--groups', required=True, help='2 column file with the sample name (1st col) and group'
                                                        ' it belongs (2nd col). Order of the samples in the file'
                                                        ' should follow the same as provided to "--b1" when running'
                                                        ' rMATS.')
    parser.add_argument('--NAs_allowed', type=int, default=0, help="Number of samples per event allowed to have"
                                                                   " NA values. Default: 0")
    parser.add_argument('--paired', action='store_true', help="Whether samples in groups file are paired.")

    rmats = parser.add_argument_group('rMATS related arguments')
    rmats.add_argument('--files', help='Events files produced by rMATS ([ev_type].MATS.[JC|JCEC].txt, one per line.'
                                       ' These files should represent a rMATS run where no statistical comparison was'
                                       ' performed ("--statoff" and "--b2" not set')

    vasttools = parser.add_argument_group('vast-tools related arguments')
    vasttools.add_argument('--tidy', help='Inclusion table produced by vast-tools tidy.')
    args = parser.parse_args()

    groups_dict = {}
    with open(args.groups) as f:
        for line in f:
            l = line.rstrip().split("\t")
            groups_dict[l[0]] = l[1]
    f.close()

    if args.tool == "rmats":
        assert args.files is not None, "When producing PCA from rMATS, please set its group arguments accordingly."
        process_rMATS(args.files, groups_dict, args.NAs_allowed, args.paired)

    elif args.tool == "vastools":
        assert args.tidy is not None, "When producing PCA from vast-tools, please set its group arguments accordingly."
        process_vasttools(args.tidy, groups_dict, args.NAs_allowed, args.paired)

if __name__ == "__main__":
    main()
