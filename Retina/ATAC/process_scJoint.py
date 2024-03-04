import pandas as pd
import numpy as np
import seaborn as sns
import umap

rna_embeddings = np.loadtxt('./output/Retina_rna_embeddings.txt')
atac_embeddings = np.loadtxt('./output/Retina_atac_embeddings.txt')
print(rna_embeddings.shape)
print(atac_embeddings.shape)
embeddings =  np.concatenate((rna_embeddings, atac_embeddings))
print(embeddings.shape)

reducer = umap.UMAP()
umap_embedding = reducer.fit_transform(embeddings)
umap_embedding.shape
df = pd.DataFrame()
df['UMAP1'] = umap_embedding[:,0]
df['UMAP2'] = umap_embedding[:,1]

rna_labels = np.loadtxt('./data/Retina_cellType_rna.txt')
atac_predictions = np.loadtxt('./output/Retina_atac_knn_predictions.txt')
labels =  np.concatenate((rna_labels, atac_predictions))
label_to_idx = pd.read_csv('./data/label_to_idx.txt', sep = '\t', header = None)
label_to_idx.shape
label_dic = []
for i in range(label_to_idx.shape[0]):
  label_dic = np.append(label_dic, label_to_idx[0][i][:-2])

data_label = np.array(["scRNA-seq", "scATAC-seq"])
df['data'] = np.repeat(data_label, [rna_embeddings.shape[0], atac_embeddings.shape[0]], axis=0)
df['predicted'] = label_dic[labels.astype(int)]

df.to_csv('./output/Retina_df.txt', index=None, sep='\t')

