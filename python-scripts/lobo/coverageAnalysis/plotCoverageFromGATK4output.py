import argparse
import os
import sys
import logging
import glob
import subprocess
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr
from collections import defaultdict, OrderedDict
import numpy as np
import math
def wccount(filename):
    out = subprocess.Popen(['wc', '-l', filename],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT
                         ).communicate()[0]
    return int(out.partition(b' ')[0])

def countOcurrences(filename,char):
    with open(filename) as f:
        return sum(char in line for line in f)

#def plotWGShisto(filename):
def processTargetedExperiment(allMetrics,perTargetMetrics):
    logging.info("Starting targeted experiment analysis.")
    if len(list(glob.glob('*HS_metrics.txt'))) > 0: #if individual metrics files exist
       logging.info("We detected individual sample metrics files. Will check the low quality base pair distribution in the data")
       processIndividualTargetSamples()

    targetedAllMetricsProcess(allMetrics)
    processPerTargetCoverage(perTargetMetrics)
def processIndividualTargetSamples():
    dict_bq={}
    nleft=len(list(glob.glob('*HS_metrics.txt')))
    for indiv_metrics in list(glob.glob('*HS_metrics.txt')):
        sname=indiv_metrics.split("_")[0]
        #logging.info("{} file! {} left..".format(indiv_metrics, nleft))
        nleft -= 1
        with open(indiv_metrics, "r") as f:
            baseQualities=[]
            for line in f:
                if "coverage_or_base_quality" in line:
                    for i, line in enumerate(f):
                        if i < 42:
                            baseQualities.append(int(line.rstrip().split("\t")[2]))
                        else:
                            if sname in dict_bq.keys():
                                logging.error("Repeated {} sample in directory.".format(indiv_metrics))
                                exit(1)
                            else:
                                dict_bq[sname] = baseQualities
                                break
        f.close()
    boxplot_dic={}
    for sample, baseQualities in dict_bq.items():
        totalbp=sum(baseQualities)
        fraction_5=sum(baseQualities[0:5])/totalbp*100
        fraction_10=sum(baseQualities[0:10])/totalbp*100
        fraction_15=sum(baseQualities[0:15])/totalbp*100
        fraction_20=sum(baseQualities[0:20])/totalbp*100
        boxplot_dic[sample] = [fraction_5,fraction_10,fraction_15,fraction_20]

    #draw boxplots
    sns.set_style("whitegrid")
    df=pd.DataFrame.from_dict(boxplot_dic)
    df=df.rename({0: "qual_5", 1: "qual_10", 2: "qual_15", 3: "qual_20"}, axis='index').transpose()
    sns.boxplot(data=df)
    plt.xlabel('Base qualities below given values')
    plt.ylabel('Percentage of base pairs (%)')
    plt.savefig("output/lowerQualityBPDistribution.pdf", format='pdf')
    plt.close()

def targetedAllMetricsProcess(infile):
    logging.info("Reading overall metrics file..")
    correl_on_off_baits=defaultdict()
    mean_median_cov=defaultdict()
    coverage_fractions=defaultdict()
    with open(infile,"r") as f:
        for line in f:
            if line.startswith("sample") and "bait_region" in line:
                header=line.rstrip().split("\t")
            else:
                (sample, bait_region, target_region, nreads, fractionBP_filtered_on_Near_baits, fractionBP_filtered_away_baits,
                 fractionBP_filt_offTarget, mean_cov_targets,median_cov_targets,fraction_targetsNoCov, fraction_tgt_1X,
                 fraction_tgt_2X, fraction_tgt_10X, fraction_tgt_20X, fraction_tgt_30X,fraction_tgt_40X,fraction_tgt_50X,fraction_100X_cov) = line.rstrip().split()
                correl_on_off_baits[sample] = [int(nreads),float(fractionBP_filtered_on_Near_baits)]
                mean_median_cov[sample] = [float(mean_cov_targets),float(median_cov_targets)]
                l=[float(i) for i in [fraction_tgt_1X, fraction_tgt_2X, fraction_tgt_10X, fraction_tgt_20X, fraction_tgt_30X, fraction_tgt_40X, fraction_tgt_50X,fraction_100X_cov]]
                coverage_fractions[sample] = l
                if float(fraction_targetsNoCov) > 0:
                    logging.info("WARNING, sample {} has some targets ({}) that didn't reach coverage=1 over any base (fully uncaptured)".format(sample,fraction_targetsNoCov))
        #violin_plot_fractionBP_on_near_baits
        logging.info("Plotting violin plot of the bait efficiency")
        df = pd.DataFrame.from_dict(correl_on_off_baits,orient="index").rename({0: "number_of_reads", 1: "fraction_bp_on_near_baits"}, axis ='columns')
        df[['fraction_bp_on_near_baits']].apply(pd.to_numeric)
        badsamples = df['fraction_bp_on_near_baits'] < 0.5
        if len(df[badsamples].index) > 0:
            logging.info("WARNING. {} samples have a high fraction (> 50%) of their aligned bases outside the baits intervals (including 250bp outside "
                         "their boundaries)".format(len(df[badsamples].index)))

            #outf.write("#List of samples with a small fraction of aligned based within the baits intervals.")
            df_bad=df[badsamples].sort_values('fraction_bp_on_near_baits', ascending=False)
            df_bad.to_csv("output/bad_wet-lab_efficiency.txt", sep='\t', columns=list(df_bad),encoding='utf-8')

        sns.set()
        sns.violinplot(data=df[df.columns[1]])
        plt.ylim(0,1)
        plt.ylabel("Fractions")
        plt.xlabel("Aligned base-pairs within baits intervals.")
        plt.axes().get_xaxis().set_ticks([])
        plt.title("Distribution of the targeted experiment efficiency across samples")
        plt.savefig("output/targetedExperimentEfficiencyDistribution.pdf", format='pdf')
        plt.close()


        ##scatterPlot - correlation between the number of reads and the bait efficiency
        logging.info("Plotting correlation between sequencing throughput and bait efficiency")
        pcc=pearsonr(df["number_of_reads"],df["fraction_bp_on_near_baits"])
        logging.info("\tPearson correlation coefficient: {}".format(pcc[0]))
        sns.regplot(x=df["number_of_reads"], y=df["fraction_bp_on_near_baits"], scatter=True,color="slategray")
        plt.title("Correlation between the number of reads and bait efficiency")
        plt.savefig("output/throughputVsBaitEfficiency.pdf",format='pdf')
        plt.close()

        #mean distribution
        logging.info("Analysing targets coverage")
        df_mean = pd.DataFrame.from_dict(mean_median_cov, orient="index").rename({0: "Mean", 1: "Median"},axis='columns')
        lowcoverage = df_mean['Median'] < 30
        if len(df_mean[lowcoverage]) > 0:
            logging.info("WARNING. {} samples were flagged with the low coverage flag because their median coverage depth value across the targets was lower than 30".format(len(df_mean[lowcoverage].index)))
            with open("output/low_averageCoverageDepth_on_target.txt", "w") as outf:
                df_lowcov=df_mean[lowcoverage].sort_values('Median', ascending=False)
                df_lowcov.to_csv("output/low_averageCoverageDepth_on_target.txt", sep='\t', columns=list(df_lowcov),encoding='utf-8')
            outf.close()
        sns.set_style("whitegrid", {'axes.grid': False})
        sns.distplot(df_mean["Median"],bins=10)
        sns.distplot(df_mean["Mean"],bins=10)
        plt.xlim(0,df_mean.iloc[:,0].max() + 50)
        #plt.title("Distribution of the mean and median across multiple samples.")

        plt.xlabel("Read depth")
        plt.ylabel("Frequency")
        plt.legend(['Median', 'Mean'])
        plt.savefig("output/coverageDepthDistributionOnTargets.pdf",format='pdf')
        plt.close()


        #fraction of targets with X coverage
        #draw boxplots
        sns.set_style("whitegrid")
        df_fractionXCoverage = pd.DataFrame.from_dict(coverage_fractions,orient="index").rename({0:"1X", 1:"2X", 2:"10X", 3:"20X", 4:"30X", 5:"40X", 6:"50X", 7:"100x"},axis='columns')
        lowgenomeCoverageOverall = df_fractionXCoverage['10X'] < 0.5
        if len(df_fractionXCoverage[lowgenomeCoverageOverall]) > 0:
            logging.info("WARNING. {} samples have more than 50% of the target regions with coverage lower than 10X.".format(len(df_fractionXCoverage[lowgenomeCoverageOverall].index)))
            with open("output/low_targetCoverageBreadth.txt", "w") as outf:
                df_lowtargetcov=df_fractionXCoverage[lowgenomeCoverageOverall].sort_values('2X', ascending=False)
                df_lowtargetcov.to_csv("output/low_targetCoverageBreadth.txt", sep='\t', columns=list(df_lowtargetcov),encoding='utf-8')
            outf.close()

        sns.boxplot(data=df_fractionXCoverage)
        plt.xlabel('Coverage')
        plt.ylabel('Fraction of target bases')
        plt.savefig("output/targetCoverageBreadthDistribution_boxplot.pdf", format='pdf')
        plt.close()

        i=1

        if len(df_fractionXCoverage.index) > 10:
            leg_value=False
        else:
            leg_value=True

        df_fractionXCoverage.T.plot(legend=leg_value,xlim=(0,7),ylim=(0,1.0),grid=True)
        x=list(df_fractionXCoverage.T.index)
        xi = [i for i in range(0, len(x))]
        plt.xticks(xi,x)
        plt.ylabel("Fraction of the target regions covered")
        plt.xlabel("Coverage depth")
        plt.savefig("output/targetCoverageBreadthDistribution_linePlot.pdf",format='pdf')
        plt.close()
def processPerTargetCoverage(perTargetMetrics):
    logging.info("Processing perTarget coverage file..")
    perTarget_dict=OrderedDict()
    sample=""
    unique_feature=set()
    perSample_dict,perSubTargetAverageCov=defaultdict(list),defaultdict(list)
    with open(perTargetMetrics, "r") as infile:
        for line in infile:
            if line.startswith("##"):
                if len(perSubTargetAverageCov) > 0:
                    perSample_dict[sample] = perSubTargetAverageCov
                sample=line[2:].rstrip()
                perSubTargetAverageCov=defaultdict(list)
            elif not line.startswith("chrom"):
                (chrom,start,end,featureLength,name,gc,mean_cov,normalized_cov,min_normalized_cov,max_normalized_cov,min_cov,max_cov,pct_0x,readCount) = line.rstrip().split()
                feature="{}:{}".format(name,start)
                unique_feature.add(name)
                perSubTargetAverageCov[name].append((int(featureLength),float(mean_cov), float(normalized_cov)))
                if feature in perTarget_dict.keys():
                    perTarget_dict[feature].append((sample,int(featureLength),float(gc),float(mean_cov),float(normalized_cov),float(min_normalized_cov),float(max_normalized_cov),int(min_cov),int(max_cov),float(pct_0x),int(readCount)))
                else:
                    perTarget_dict[feature] = [(sample,int(featureLength),float(gc),float(mean_cov),float(normalized_cov),float(min_normalized_cov),float(max_normalized_cov),int(min_cov),int(max_cov),float(pct_0x),int(readCount))]

        number_of_samples = set([len(v) for k,v in perTarget_dict.items()])
        if len(number_of_samples) != 1:
            logging.info("Some samples may not have results displayed for all the target features: {}".format(number_of_samples))
        else:
            logging.info("{} samples under analysis.".format(number_of_samples.pop()))



        # #sub dicts for specific plots
        # #perFeature coverage vs gc content scatterplot
        df_corr_cov_gc= pd.DataFrame()
        # #perFeature coverage vs feature length
        df_corr_cov_length=pd.DataFrame()

        #per sample and perFeature_normalizedCov
        first_val=list(perTarget_dict.values())[0]
        cols=[sample[0] for sample in first_val]
        dfperFeatureAndPerSample_norm_cov=pd.DataFrame(index=perTarget_dict.keys(), columns=set([sample[0] for sample in first_val]))

        #features with no coverage
        df_featuresWith0X=pd.DataFrame()
        for feature, values in perTarget_dict.items():
            perFeature_avg_cov=sum([e[4] for e in values])/len(values)
            perFeature_avg_0cov=sum([e[9] for e in values])/len(values)
            df2 = pd.DataFrame({"gc_content" : [values[0][2]], "mean_normalized_coverage" : [perFeature_avg_cov]}, index=[feature])
            df3 = pd.DataFrame({"feature_length" : [values[0][1]], "mean_normalized_coverage": [perFeature_avg_cov]}, index=[feature])

            for sample in values:
                dfperFeatureAndPerSample_norm_cov.at[feature,sample[0]] = sample[4]

            df5 = pd.DataFrame({"mean_fraction_cov0" : [perFeature_avg_0cov]}, index=[feature])
            df_corr_cov_gc = df_corr_cov_gc.append(df2)
            df_corr_cov_length = df_corr_cov_length.append(df3)
            df_featuresWith0X=df_featuresWith0X.append(df5)

        p10=df_corr_cov_gc['mean_normalized_coverage'].quantile(0.1)
        logging.info("Analysing mean coverage per target")
        logging.info("Mean normalized coverage value for the 10th percentile (10% of the features have mean coverage below this value, which translates for targets with lower coverage): {}".format(p10))
        logging.info("Writing targets with lowest coverage to 'targetsWithLowestMeanCoverage.txt' file")
        df_corr_cov_gc[df_corr_cov_gc.mean_normalized_coverage < df_corr_cov_gc['mean_normalized_coverage'].quantile(0.1)].to_csv("output/targetsWithLowestMeanCoverage.txt",sep="\t", columns=list(df_corr_cov_gc))

        ###scatterplots
        logging.info("Plotting per target correlation between gc content and mean coverage")
        pcc=pearsonr(df_corr_cov_gc["gc_content"],df_corr_cov_gc["mean_normalized_coverage"])
        logging.info("\tPearson correlation coefficient: {}".format(pcc[0]))
        sns.regplot(x=df_corr_cov_gc["gc_content"], y=df_corr_cov_gc["mean_normalized_coverage"], scatter=True,color="slategray")
        plt.title("Correlation between gc content and the mean coverage")
        plt.savefig("output/perTarget_gcVscoverage.pdf",format='pdf')
        plt.ylim(0,3)
        plt.close()

        logging.info("Plotting per target correlation between feature length and mean coverage")
        pcc=pearsonr(df_corr_cov_length["feature_length"],df_corr_cov_length["mean_normalized_coverage"])
        logging.info("\tPearson correlation coefficient: {}".format(pcc[0]))
        sns.regplot(x=df_corr_cov_length["feature_length"], y=df_corr_cov_length["mean_normalized_coverage"], scatter=True,color="green")
        plt.title("Correlation between feature length and mean coverage")
        plt.savefig("output/perTarget_featureLengthVscoverage.pdf",format='pdf')
        plt.close()
        ##heatmap
        #lt.pcolor(dfperFeatureAndPerSample_norm_cov.apply(pd.to_numeric,errors='coerce'))
        #sns.heatmap(# .apply(pd.to_numeric, errors='coerce'))
        #lt.show()
        #targets not covered
        logging.info("Analysing fractions of targets with 0 coverage.")
        logging.info("Plotting per target average of regions with no coverage")
        bins_count = df_featuresWith0X.groupby(pd.cut(df_featuresWith0X.mean_fraction_cov0, right=True,bins=[ -0.01,0,0.1,0.2,0.5,0.75,1], include_lowest=False)).count()
        bins_count.plot(kind='bar',legend=False,color=sns.color_palette("GnBu_d"))
        plt.xticks([0,1,2,3,4,5,6], ["100%","> 90% ", "> 80%", ">50%", ">25%", ">0%" ])
        plt.tick_params(axis='x', rotation=0)
        plt.ylabel("Number of features")
        plt.xlabel("Average coverage breadth")
        plt.savefig("output/perTarget_averageCoverageBreadth.pdf",format='pdf')
        plt.close()
        df_featuresWith0X[df_featuresWith0X["mean_fraction_cov0"] > 0.1].sort_values(by="mean_fraction_cov0",ascending=False).to_csv("output/perTarget_withSignificantFractionNotCovered.txt",sep="\t",columns=list(df_featuresWith0X))

        logging.info("Testing for coding/non coding coverage analysis.")
        df_codingNonCodingCov = pd.DataFrame(index=perSample_dict.keys(),columns=unique_feature)
        with open("output/perTargetInfoAbout0XCoverageSubfeatures.txt", "w") as outf:
            for sample, feature in perSample_dict.items():
                outf.write("##{}\n".format(sample))
                for feat_name, subfeature in feature.items():
                    cov_0X=[i[1] for i in subfeature if i[1] == 0.0]
                    if len(cov_0X) > 0:
                        outf.write("{} feature has {} subfeatures, where {} if them has no coverage at all.\n".format(feat_name,len(subfeature),len(cov_0X)))
                    avg_cov=np.average([i[1] for i in subfeature],weights=[j[0] for j in subfeature])
                    df_codingNonCodingCov.at[sample, feat_name] = avg_cov
        outf.close()

        df_codingNonCodingCov.to_csv("output/exon.csv", sep="\t")
        #df.columns = [i if "_" not in i else i + "=" + str(newElements[int(i[-1]) - 1]) for i in df.columns]
        df_codingNonCodingCov.columns = [i.split("_")[1] for i in df_codingNonCodingCov.columns if "exon" or "intron" in i ]
        #print(df_codingNonCodingCov)

        plt.figure(figsize=(9,8))
        sns.boxplot(data=df_codingNonCodingCov)
        plt.xticks(rotation='vertical')
        #plt.ylim(0,1500)
        plt.ylabel("Average read depth")
        plt.savefig("output/perTarget_depthOfCoverage.pdf", format="pdf")
def main():
    parser = argparse.ArgumentParser(description='Script to plot useful data from a genome coverage analysis performed with GAT4K4. Matplotlib is required')
    parser.add_argument(dest='path',help='Path to the directory which holds all the output files.')
    parser.add_argument(dest='analysis_type', metavar='analysis_type', choices=("WGS", "targeted"), help="Type of the analysis performed. Choices: [WGS,targeted]")
    args = parser.parse_args()
    try:
        os.chdir(args.path)
        if not os.path.exists("output"):
            os.makedirs("output")
        if args.analysis_type == "WGS":
            wgs_all="final_colllectWgsMetrics_all.txt"
            if not os.path.exists(wgs_all):
                logging.error("Error: '{}' file doesn't exist in {} directory, which is required when WGS analysis is set.".format(wgs_all, args.path))
                sys.exit(1)
            elif len(glob.glob('*.histo')) == 0:
                logging.error("Error: No histogram files '*.histo' in {} directory. They are required in WGS experiments.".format(args.path))
                sys.exit(1)
            elif len(glob.glob('*.histo')) != (wccount(wgs_all) - 1):
                logging.error("Error. Number of '.histo' files must be the same as the number of samples described in the {} file.".format(wgs_all))
                sys.exit(1)
            #else:
            #    for f in glob.glob('*.histo'):

        elif args.analysis_type == "targeted":
            targeted_all = "final_collectHSmetrics_all.txt"
            perTarget_all = "final_perTargetCoverage.txt"
            if not os.path.exists(targeted_all) or not os.path.exists(perTarget_all):
                logging.error("Error: '{}' and '{}' must exist in {} directory. They are required when targeted analysis is set.".format(targeted_all, perTarget_all,args.path))
                sys.exit(1)
            elif countOcurrences(perTarget_all, "##") != (wccount(targeted_all) - 1):
                logging.error("Error. Number of samples differs in {} and {} files. ({};{})".format(targeted_all, perTarget_all,wccount(targeted_all)-1,countOcurrences(perTarget_all,"##")))
                sys.exit(1)
            else:
                processTargetedExperiment(targeted_all, perTarget_all)

    except OSError:
        logging.error("Error: {} is not a directory!".format(args.path))

if __name__ == "__main__":
    main()
