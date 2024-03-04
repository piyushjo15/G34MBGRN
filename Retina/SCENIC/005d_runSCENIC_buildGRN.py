
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

projDir = "Retina/SCENIC/"
tmpDir = 'TMP'

#load Scenicplus object
import pickle
infile = open(os.path.join(projDir, 'scenicplus',Sample,'scplus_obj.pkl'), 'rb')
scplus_obj = pickle.load(infile)
infile.close()

############# Build eGRNs ##############

# Load functions
from scenicplus.grn_builder.gsea_approach import build_grn


build_grn(scplus_obj,
         min_target_genes = 10,
         adj_pval_thr = 1,
         min_regions_per_gene = 0,
         quantiles = (0.85, 0.90, 0.95),
         top_n_regionTogenes_per_gene = (5, 10, 15),
         top_n_regionTogenes_per_region = (),
         binarize_using_basc = True,
         rho_dichotomize_tf2g = True,
         rho_dichotomize_r2g = True,
         rho_dichotomize_eregulon = True,
         rho_threshold = 0.05,
         keep_extended_motif_annot = True,
         merge_eRegulons = True,
         order_regions_to_genes_by = 'importance',
         order_TFs_to_genes_by = 'importance',
         key_added = 'eRegulons_importance',
         cistromes_key = 'Unfiltered',
         disable_tqdm = False, #If running in notebook, set to True
         ray_n_cpu = 4,
         _temp_dir = tmpDir)

print("Finished Build GRN")
from scenicplus.utils import format_egrns
format_egrns(scplus_obj, eregulons_key = 'eRegulons_importance', TF2G_key = 'TF2G_adj', key_added = 'eRegulon_metadata')
print("check unstructured data in scplus object after adding eRegulon metadata")
scplus_obj.uns.keys()

directory_path = os.path.join(projDir,'scenicplus/',Sample,'output/')

# Check if the directory already exists
if not os.path.exists(directory_path):
    os.makedirs(directory_path)


directory_path = os.path.join(projDir,'scenicplus/',Sample,'output/metadata/')

# Check if the directory already exists
if not os.path.exists(directory_path):
    os.makedirs(directory_path)

############# save metadata and  AUC matrix #####################

metadata = scplus_obj.uns['eRegulon_metadata']
metadata.to_csv('Retina/SCENIC/scenicplus/' + Sample + '/output/' +Sample +'_eRegulon_metadata.csv',
                float_format='%.3f',index=False, sep=';')


print("Saved eRegulon metadata for sample:" + Sample)

# Save
#to access the eGRNs
import dill
with open(os.path.join(projDir, 'scenicplus/',Sample,'scplus_obj.pkl'), 'wb') as f:
  dill.dump(scplus_obj, f)

print("saved after bulding eGRNs for Sample:" +Sample)


import dill
with open(os.path.join(projDir, 'scenicplus/',Sample,'scplus_obj_backup.pkl'), 'wb') as f:
  dill.dump(scplus_obj, f)
print("Saved scplus_obj backup!")

