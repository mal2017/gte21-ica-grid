from random import randint
from time import sleep

# running into ftputil errors, one of the below libraries must check with an ftp site
# for updates. TO deal with this I have the script sleep a random interval
# after each module loads. This seems to massively improve the issue.
# Ideally there would be a slurm flag that makes sure jobs aren't submitted
# all at once, but I can't find anything like it.
sleep(randint(1,60))

import pandas as pd
sleep(randint(1,60))
import numpy as np
sleep(randint(1,60))
import random
sleep(randint(1,60))
from sklearn.decomposition import FastICA

seed = snakemake.params['random_seed']
random.seed(seed)
np.random.seed(seed)

X = pd.read_csv(snakemake.input[0])
X.set_index('index',inplace=True)
X = X.T
ica = FastICA(n_components=int(snakemake.params['comps']), random_state=seed, max_iter=5000)
source = ica.fit_transform(X)

mixing = pd.DataFrame(ica.mixing_, X.columns)
res = pd.DataFrame(source, X.index)

mixing.to_csv(snakemake.output["mixing"])
res.to_csv(snakemake.output["source"])
