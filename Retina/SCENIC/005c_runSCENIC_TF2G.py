
###### add TF2G Adjancencies #########
import scenicplus
scenicplus.__version__

import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')

#define dir
Sample = sys.argv[1]
print("Add TF2G adjancecies to: " + Sample)

projDir = "Retina/SCENIC/"
tmpDir = 'TMP/'

#load Scenicplus object
import pickle
infile = open(os.path.join(projDir, 'scenicplus',Sample,'scplus_obj.pkl'), 'rb')
scplus_obj = pickle.load(infile)
infile.close()

#load TSV file

from scenicplus.TF_to_gene import *
load_TF2G_adj_from_file(scplus_obj,
                        f_adj = projDir+Sample+'/cor'+Sample+'_cleaned.tsv',
                        inplace = True,
                        key= 'TF2G_adj')

print('Checking if TF2G has been added to scplus_obj')
print(scplus_obj.uns.keys())

# Save
# Save
import pickle
pickle.dump(scplus_obj,
            open(os.path.join(projDir, 'scenicplus',Sample,'scplus_obj.pkl'), 'wb'))
  
print("Saved scenicplus object after adding TF2G adjacencies for "+Sample)


