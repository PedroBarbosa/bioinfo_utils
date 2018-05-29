import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from collections import defaultdict
import numpy as np
df_exon=pd.read_csv(sys.argv[1],sep="\t")
for i in df_exon.columns:
    print(i.split("_")[1])
genes=[i.split("_")[1] for i in df_exon.columns]
print(genes)
df_exon.rename(columns={genes})
df_intron=pd.read_csv(sys.argv[2],sep="\t")
#print(df_intron.head())
df_merge=pd.concat([df_exon,df_intron])
df=df_merge.stack()
print(df.head())

genes=set([i.split("_")[1] for i in df_exon.columns if "exon" or "intron" in i ])
a=df_exon.to_dict()
b=df_intron.to_dict()
final_dict_e=defaultdict(list)
final_dict_i=defaultdict(list)
for g in genes:
    for feat in a.keys():
        if g == feat.split("_")[1]:
            final_dict_e[g].append(a[feat].values())
            final_dict_i[g].append(b[feat.replace("exon","intron")].values())

df_e=pd.DataFrame.from_dict(final_dict_e)
df_e.columns=list(genes)
df_i=pd.DataFrame.from_dict(final_dict_i)
df_i.columns=list(genes)

df_fin=pd.concat([df_i,df_e])


#print(df_fin["MYH7"])
#df_merge.index[-1]
#print(df_merge)
#plt.figure(figsize=(9, 8))
#sns.boxplot(data=df_fin)
#sns.despine(offset=10, trim=True)
#plt.xticks(rotation='vertical')
#plt.show()
#df_merge["gene"] = [i.split("_")[1] for i in df_merge.columns if "exon" or "intron" in i ][0]
#print(df_merge["gene"].head())
