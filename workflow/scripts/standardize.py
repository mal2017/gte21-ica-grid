import pandas as pd
import numpy as np

def standardize(X, maxval):
    Z = (X - X.mean())/(X.std())
    if maxval:
      return Z.clip(-1*maxval,maxval)
    else:
      return Z

expression = pd.read_csv(snakemake.input['expr'])
genes = pd.read_csv(snakemake.input['genes']).genes
expression.set_index('index', inplace=True)
expression = expression.loc[:,genes]
#expression = expression.T
expression = standardize(expression, snakemake.params["maxval"]).T
expression.reset_index().to_csv(snakemake.output[0], index=False)
