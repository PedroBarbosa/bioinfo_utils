import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from collections import defaultdict
import numpy as np



df = pd.DataFrame(columns=["coverage", "gene","type"])
for f in range(1,3):
    if f == 1:
        type="exon"
    else:
        type="intron"
    feat_dict = {}
    feat_name = ""
    aux = []
    print("Reading file {}".format(f))
    with open(sys.argv[f],'r') as infile:
        for line in infile:
            if not line.startswith("chrom") and not line.startswith("##"):
                if line.split()[4] != feat_name:
                    if not feat_name == "":

                        weighted_avg = np.average([x[0] for x in aux], weights=[x[1] for x in aux])
                        if feat_name.split("_")[1] in feat_dict.keys():
                            feat_dict[feat_name.split("_")[1]].append(weighted_avg)
                        else:
                            feat_dict[feat_name.split("_")[1]] = [weighted_avg]

                    feat_name=line.split()[4]
                    aux=[(float(line.split("\t")[6]),int(line.split("\t")[3]))]

                else:
                    aux.append((float(line.split("\t")[6]),int(line.split("\t")[3])))

        if not feat_name == "":
            weighted_avg = np.average([x[0] for x in aux], weights=[x[1] for x in aux])
            if feat_name.split("_")[1] in feat_dict.keys():
                feat_dict[feat_name.split("_")[1]].append(weighted_avg)
            else:
                feat_dict[feat_name.split("_")[1]] = weighted_avg


    infile.close()
    print("Appending to dataframe..")
    for k,v in feat_dict.items():
        for val in v:
            df = df.append({
                "coverage": val,
                "gene": k,
                "type": type
            }, ignore_index=True)


df.to_csv("df.ready.csv")
df=pd.read_csv("df.ready.csv")
print(df.shape)
#df_=pd.concat([df,df2])



plt.figure(figsize=(15,8))

sns.set(style="ticks")
sns.boxplot(x="gene",y="coverage", data=df,hue="type",saturation=0.8,width=0.8,palette="PRGn")
#sns.despine(offset=10, trim=True)
plt.xticks(rotation='vertical')
plt.show()
#plt.savefig("perTarget_depthOfCoverage.pdf", format="pdf")