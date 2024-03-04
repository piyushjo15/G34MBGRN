
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
import numpy as np

_stderr = sys.stderr
null = open(os.devnull,'wb')


#define dir and Sample
import sys
Sample = sys.argv[1]
print("Start Analysis for Sample: " + Sample)
projDir = "Retina/SCENIC/"
tmpDir = 'TMP/'


#Create Scenic object #############
# Load functions
from scenicplus.scenicplus_class import SCENICPLUS, create_SCENICPLUS_object
from scenicplus.preprocessing.filtering import *
  
  #load scRNA and scATAC data
  # Load data
  
adata = sc.read_h5ad(os.path.join(projDir,Sample,'adata.h5ad'))
cistopic_obj = dill.load(open(os.path.join(projDir, Sample,'cistopic_obj.pkl'), 'rb'))
menr = dill.load(open(os.path.join(projDir, Sample,'pycistarget','motifs/menr.pkl'),'rb'))

#Create SCENIC+ object - for Non-Multiome data
from scenicplus.scenicplus_class import create_SCENICPLUS_object
scplus_obj = create_SCENICPLUS_object(
        GEX_anndata = adata,
        cisTopic_obj = cistopic_obj,
        menr = menr,
        multi_ome_mode = False,
        key_to_group_by = 'predictedGroup',
        nr_cells_per_metacells = 5)


#scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense()) ##already nparray
scplus_obj


print(f"The cell lines for which we have scRNA-seq data are:\t{', '.join(set(adata.obs['predictedGroup']) - set(['-']))}")
print(f"The cell lines for which we have scATAC-seq data are:\t{', '.join(set(cistopic_obj.cell_data['predictedGroup']))}")
print(f"The cell lines for which we have both:\t{', '.join(set(cistopic_obj.cell_data['predictedGroup']) & set(adata.obs['predictedGroup']))}")


# Save
if not os.path.exists(os.path.join(projDir, 'scenicplus',Sample)):
  os.makedirs(os.path.join(projDir, 'scenicplus',Sample))


import pickle
pickle.dump(scplus_obj,
            open(os.path.join(projDir, 'scenicplus',Sample,'scplus_obj.pkl'), 'wb'))

print("Saved scplus object")

#Check with which Biomart host the gene names used in analysis match best
ensembl_version_dict = {'105': 'http://www.ensembl.org',
  '104': 'http://may2021.archive.ensembl.org/',
  '103': 'http://feb2021.archive.ensembl.org/',
  '102': 'http://nov2020.archive.ensembl.org/',
  '101': 'http://aug2020.archive.ensembl.org/',
  '100': 'http://apr2020.archive.ensembl.org/',
  '99': 'http://jan2020.archive.ensembl.org/',
  '98': 'http://sep2019.archive.ensembl.org/',
  '97': 'http://jul2019.archive.ensembl.org/',
  '96': 'http://apr2019.archive.ensembl.org/',
  '95': 'http://jan2019.archive.ensembl.org/',
  '94': 'http://oct2018.archive.ensembl.org/',
  '93': 'http://jul2018.archive.ensembl.org/',
  '92': 'http://apr2018.archive.ensembl.org/',
  '91': 'http://dec2017.archive.ensembl.org/',
  '90': 'http://aug2017.archive.ensembl.org/',
  '89': 'http://may2017.archive.ensembl.org/',
  '88': 'http://mar2017.archive.ensembl.org/',
  '87': 'http://dec2016.archive.ensembl.org/',
  '86': 'http://oct2016.archive.ensembl.org/',
  '80': 'http://may2015.archive.ensembl.org/',
  '77': 'http://oct2014.archive.ensembl.org/',
  '75': 'http://feb2014.archive.ensembl.org/',
  '54': 'http://may2009.archive.ensembl.org/'}

import pybiomart as pbm
def test_ensembl_host(scplus_obj, host, species):
  dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
  annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
  annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
  annot['Chromosome'] = annot['Chromosome'].astype('str')
  filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
  annot = annot[~filter]
  annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
  gene_names_release = set(annot['Gene'].tolist())
  ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
  print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
  return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
  print(f'host: {version}')
  try:
    n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'hsapiens')
  except:
    print('Host not reachable')

v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")

