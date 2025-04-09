import os
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RNAY8= pd.read_csv("/~/PFCdata/RNAY8/RNAY8_raw.txt", sep='\t', header=0, index_col=0).T
scrub = scr.Scrublet(RNAY8, expected_doublet_rate=0.06,sim_doublet_ratio=20)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, 
                                                          min_cells=0, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
predicted_doublets=scrub.call_doublets(threshold = 0.15)
os.chdir('/~/PFCdata/RNAY8/')
np.savetxt('RNAY8_doublet_position.csv', predicted_doublets, delimiter="\t")

