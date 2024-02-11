## This script runs topic modeling on cisTopic object

import os
import pandas
import numpy
import warnings
import pycisTopic
import sys

#Define dir
Sample = sys.argv[1]
print("Start Analysis for Sample: " + Sample)

projDir = os.path.join("SCENIC/scATAC/",Sample)

tmpDir = 'TMP/'

import pickle
cistopic_obj = pickle.load(open(os.path.join(projDir, 'cistopic_obj.pkl'), 'rb'))

from pycisTopic.cistopic_class import *
models=run_cgs_models(cistopic_obj,
                        n_topics=[50],
                        n_cpu=4,
                        n_iter=500,
                        random_state=555,
                        alpha=50,
                        alpha_by_topic=True,
                        eta=0.1,
                        eta_by_topic=False,
                        save_path=None,
                        _temp_dir = tmpDir)

if not os.path.exists(os.path.join(projDir, 'models')):
  os.makedirs(os.path.join(projDir, 'models'))

pickle.dump(models,
            open(os.path.join(projDir, 'models/ATAC_models_500_iter_LDA.pkl'), 'wb'))


print("Finished Topic Modeling and saved Models for Sample" + Sample)


exit()


