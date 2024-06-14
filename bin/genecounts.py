#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import sys
from functools import reduce

args = sys.argv[1:]
directory = ''.join(args)

filenames = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.endswith('.featurecounts')]

dfs= []
all_counts = filenames
for file in all_counts:
    f = pd.read_table((directory + "/" + file), comment="#")
    f1 = f.iloc[:,[0,6]]
    fc_colname = f1.columns[1]
    fc_colname_new = fc_colname.split(sep="/")[-1]
    fc_colname_newer = fc_colname_new.replace("_Aligned.out.bam", "")
    f2= f1.rename(columns = {fc_colname: fc_colname_newer})
    dfs.append(f2)

def merge_df(a, b):
    return pd.merge(a, b, how = "outer", on = "Geneid")

df_final = reduce(merge_df,dfs)
out = "genecounts.csv"
df_final.to_csv(out, index = False)
