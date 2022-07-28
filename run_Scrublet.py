## module add Python/3.8.6-GCCcore-10.2.0
## Detect doublets using scrublet
## https://github.com/swolock/scrublet/blob/master/examples/scrublet_basics.ipynb
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import csv
import pandas as pd
# Python program to demonstrate
# command line arguments
import getopt, sys 
# Remove 1st argument from the
# list of command line arguments
argumentList = sys.argv[1:]
# options
options = "hi:o:"
# Long options
long_options = ["Help", "Input =", "Output ="]
try:
    # Parsing argument
    arguments, values = getopt.getopt(argumentList, options, long_options)
    # checking each argument
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-h", "--Help"):
            print ("-i input dir for filtered_feature_bc_matrix from 10x cellrangeri, must have the file matrix.mtx.gz")
             
        elif currentArgument in ("-i", "--Input"):
            inputdir=currentValue; 
        elif currentArgument in ("-o", "--Output"):
            outputname=currentValue;

except getopt.error as err:
    # output error, and return with an error code
    print (str(err))

## plotting setting
## plt.rcParams['font.family'] = 'sans-serif'
## plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42
## file dir
input_dir = inputdir;
## count matrix
counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
## gene list
#genes = np.array(scr.load_genes(input_dir + '/genes.tsv', delimiter='\t', column=1))
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
#print('Number of genes in gene list: {}'.format(len(genes)))
## defualt settings
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.076)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,min_cells=3,min_gene_variability_pctl=85,n_prin_comps=30)
## plot
scrub.plot_histogram();
plt.savefig((outputname+"_scrublet.pdf"))
## save results
df = pd.DataFrame({"predicted_doublets" : predicted_doublets, "doublet_scores" : doublet_scores})
df.to_csv((outputname+"_scrublet.csv"), index=False)

scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True);
plt.savefig((outputname+"_umap.pdf"))
