


import scenicplus
scenicplus.__version__

import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

#define dir
Sample = sys.argv[1]
print("Start Analysis for Sample: " + Sample)

projDir = os.path.join("Retina/SCENIC/")
tmpDir = 'TMP'
#load Scenicplus object
import pickle
infile = open(os.path.join(projDir, 'scenicplus',Sample,'scplus_obj.pkl'), 'rb')
scplus_obj = pickle.load(infile)
infile.close()


#Generate Cistromes
from scenicplus.cistromes import *
  
import time
start_time = time.time()
merge_cistromes(scplus_obj)
time = time.time()-start_time
print(time/60)

##### Infer Enhancer - Gene Relationships #######
#Get Search Space
biomart_host = "http://feb2021.archive.ensembl.org/" #paste correct biomart_host (defined in 004_create_scplus_obj.py script)

from scenicplus.enhancer_to_gene import get_search_space, calculate_regions_to_genes_relationships, RF_KWARGS
get_search_space(scplus_obj,
                 biomart_host = biomart_host,
                 species = 'hsapiens',
                 assembly = 'hg38',
                 upstream = [1000, 150000],
                 downstream = [1000, 150000])

#Enhancer to Gene Models using correlation = random forest or Gradient Boosting Machine

calculate_regions_to_genes_relationships(scplus_obj,
                                         ray_n_cpu = 2,
                                         _temp_dir = os.path.join(tmpDir,'ray_spill'),
                                         importance_scoring_method = 'RF',
                                         importance_scoring_kwargs = {'n_jobs': 3, 'max_features': 0.3, 'n_estimators': 1000})

# Save
import pickle
pickle.dump(scplus_obj,
            open(os.path.join(projDir, 'scenicplus',Sample,'scplus_obj.pkl'), 'wb'))
# Save backup
import pickle
pickle.dump(scplus_obj,
            open(os.path.join(projDir, 'scenicplus',Sample,'scplus_obj_backup.pkl'), 'wb'))


print("Saved scenicplus object after calculating regions-genes reslationships"+Sample)


