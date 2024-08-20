# Load required modules
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import csv
import os

os.chdir('/lustre/scratch126/tol/teams/lawniczak/users/jr35/phd/Mali2/')

## Read the donor names from the CSV file
don_nms = [line.strip() for line in open('data/processed/Pf/don_nms.txt', 'r')]

print(don_nms)

## Read in matrices and genes and put into list
counts_msc = []

for dn in don_nms:
    x = scipy.io.mmread('data/processed/Pf/' + dn + '/seu_obj/counts_msc.mtx')
    counts_msc.append(x)

# print(counts_msc[1])

genes_msc = []
for dn in don_nms:
    x = [line.strip() for line in open('data/processed/Pf/' + dn + '/seu_obj/genes_msc.txt', 'r')]
    genes_msc.append(x)

# print(genes_msc[1])

## Calculate doublet rate
dbr = [(matrix.shape[1]/1000)*0.008 for matrix in counts_msc]
print(dbr)

## Generate empty lists for storing scrublet results from loop below
# db_scores_ls = []
# prd_db_ls = []

## Run scrublet
for mat, gns, dbr in zip(counts_msc, genes_msc, dbr):
  counts_matrix = mat.T.tocsc()
  
  # Calculate number of counts and genes per cell
  nCount = np.sum(counts_matrix, axis = 1)
  nFeature = np.sum(counts_matrix > 0, axis = 1)
  
  # Load gene names
  genes = gns
  # Output shape
  print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
  print('Number of genes in gene list: {}'.format(len(genes)))
  print('Doublet rate: {}'.format(dbr))
  
  # Since the algorithm is robust to the expected doublet rate we will not search extensively for the best fit here
  scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=dbr)
  
  # Run scrublet
  doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
  min_cells=3,
  min_gene_variability_pctl=90,
  n_prin_comps=30,
  log_transform=True,
  mean_center=True,
  normalize_variance=True,
  synthetic_doublet_umi_subsampling = 1,
  verbose=False
  )

  csv.writer(open('data/processed/Pf/' + dn + '/seu_obj/doublet_scores.csv', 'w', newline='')).writerow(doublet_scores)
  csv.writer(open('data/processed/Pf/' + dn + '/seu_obj/predicted_doublets.csv', 'w', newline='')).writerow(predicted_doublets)
  
  
  # db_scores_ls.append(doublet_scores)
    
  # prd_db_ls.append(predicted_doublets)
  
