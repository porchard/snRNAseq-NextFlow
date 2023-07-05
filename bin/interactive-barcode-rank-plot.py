#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
import plotly.express as px
from scipy.io import mmread

MATRIX, OUTPUT = sys.argv[1:]

mtx = mmread(MATRIX)
mtx = pd.DataFrame.sparse.from_spmatrix(mtx)

umi_counts = mtx.sum()
umi_counts = umi_counts[umi_counts>0]

br = pd.DataFrame({'umis': umi_counts}).sort_values('umis', ascending=False)
br['rnk'] = range(1, len(br)+1)

fig = px.scatter(br, x="rnk", y="umis", log_x=True, log_y=True, labels={
                     "rnk": "Rank",
                     "umis": "UMIs"
                 })

fig.write_html(OUTPUT)