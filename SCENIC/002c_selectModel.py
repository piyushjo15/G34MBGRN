


import os
import pandas as pd
import numpy as np
import warnings
import pycisTopic
import pickle

#Define dir
import sys
Sample = sys.argv[1]
print("Start Analysis for Sample: " + Sample)


projDir = os.path.join("SCENIC/scATAC/",Sample)
tmpDir = 'TMP/'
####### MODEL SELECTION ##########
models = pickle.load(open(os.path.join(projDir, 'models/ATAC_models_500_iter_LDA.pkl'), 'rb'))
cistopic_obj = pickle.load(open(os.path.join(projDir, 'cistopic_obj.pkl'), 'rb'))

from pycisTopic.lda_models import *
model=evaluate_models(models,
                      select_model=50,
                      return_model=True,
                      metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                      plot_metrics=False)

# Add model to cisTopicObject
cistopic_obj.add_LDA_model(model)

# Save
pickle.dump(cistopic_obj,
            open(os.path.join(projDir, 'cistopic_obj.pkl'), 'wb'))


############## Visualization ######################
palette={"#FF6EB4","#42B86E","#A8D05A","#4C9A2A","#054907","#CDAF95","#66CD00","#FF4040","#FFD700","#FFFF00","#E9842C","#BBFFFF","#0000FF","#7F96FF","#00FFFF","#06C2AC"}

from pycisTopic.clust_vis import *
run_umap(cistopic_obj, target  = 'cell', scale=True)
#plot_metadata(cistopic_obj, reduction_name = 'UMAP', variables = ['Clusters_Combined'], show_legend=True,text_size=5, color_dictionary=palette)
#plot_metadata(cistopic_obj, reduction_name = 'UMAP', variables = ['celltype'], show_legend=True,text_size=5, color_dictionary=palette)

plot_topic(cistopic_obj, reduction_name = 'UMAP', num_columns = 5,save=os.path.join(projDir,'plot_Topics.pdf'))

#plot cluser
plot_metadata(cistopic_obj,
              reduction_name='UMAP',
              variables=['Clusters_Combined'], # Labels from RNA and new clusters
              target='cell',
              text_size=10,
              dot_size=5,
              figsize=(10,10),
              save= os.path.join(projDir,'dimensionality_reduction_label.pdf'))

###### Inferring Candidate enhancer region

# Next we will infer candidate enhancer regions by:
#   binarization of region-topic probabilites.
# calculation differentially accessibile regions (DARs) per cell type.


#First we will binarize the topics using the otsu method and by taking the top 3k regions per topic.
from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

#Run without imputation
#Run PycisTopic without imputation
from pycisTopic.diff_features import *
nonimputed_obj = CistopicImputedFeatures(imputed_acc=cistopic_obj.fragment_matrix,cell_names=cistopic_obj.cell_names,feature_names=cistopic_obj.region_names,project="cistopic_not_imputed")

# remove zero rows as done by the impute_accessibility function

zero_rows = (nonimputed_obj.mtx.sum(axis=1) == 0)
nonimputed_obj.mtx = nonimputed_obj.mtx[~np.array(zero_rows).flatten(), :]
nonimputed_obj.feature_names = list(np.array(nonimputed_obj.feature_names)[~np.array(zero_rows).flatten()])

#Next we will calculate DARs per cell type
normalized_imputed_acc_obj = normalize_scores(nonimputed_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, nonimputed_obj, variable='Clusters_Combined', var_features=variable_regions, split_pattern = '-')

#Save Results
if not os.path.exists(os.path.join(projDir, 'candidate_enhancers')):
  os.makedirs(os.path.join(projDir, 'candidate_enhancers'))

import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(projDir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(projDir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(projDir, 'candidate_enhancers/markers_dict.pkl'), 'wb'))

#save cistopic object
pickle.dump(cistopic_obj,
            open(os.path.join(projDir, 'cistopic_obj.pkl'), 'wb'))

print("Finished 002_selectModel for:" + Sample)

