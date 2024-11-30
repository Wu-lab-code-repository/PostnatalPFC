##########py
import os
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

sampleraw=os.listdir('/home/gpfs/home/wulab15/MacacaPFC/snRNA/postnatal/raw/')
sampleraw

for i in range(0,len(sampleraw)):
    PFCi= pd.read_csv(sampleraw[i], sep='\t', header=0, index_col=0).T
    scrub = scr.Scrublet(PFCi, expected_doublet_rate=0.06,sim_doublet_ratio=20)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, 
                                                              min_cells=0, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)
    predicted_doublets=scrub.call_doublets(threshold = 0.15)
    samplename=sampleraw[i][0:-8]
    print(samplename)
    print(predicted_doublets[predicted_doublets==False].shape)        
    os.chdir('/home/gpfs/home/wulab15/MacacaPFC/snRNA/postnatal/raw/')
    np.savetxt('%s_doublet_position.csv'%(samplename), predicted_doublets, delimiter="\t")
    pd.DataFrame(doublet_scores).to_csv('%s_doublet_score_test.csv'%(samplename),sep='\t',header=0)
    pd.DataFrame(PFCi.index).to_csv('%s_doublet_score_cell_name_test.csv'%(samplename),sep='\t',header=0)

##########py
import os
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

sampleraw=os.listdir('/home/gpfs/home/wulab15/MacacaPFC/snRNA/embryo/raw/')
sampleraw

for i in range(0,len(sampleraw)):
    PFCi= pd.read_csv(sampleraw[i], sep='\t', header=0, index_col=0).T
    scrub = scr.Scrublet(PFCi, expected_doublet_rate=0.06,sim_doublet_ratio=20)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, 
                                                              min_cells=0, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)
    predicted_doublets=scrub.call_doublets(threshold = 0.15)
    samplename=sampleraw[i][0:-8]
    print(samplename)
    print(predicted_doublets[predicted_doublets==False].shape)        
    os.chdir('/home/gpfs/home/wulab15/MacacaPFC/snRNA/embryo/raw/')
    np.savetxt('%s_doublet_position.csv'%(samplename), predicted_doublets, delimiter="\t")
    pd.DataFrame(doublet_scores).to_csv('%s_doublet_score_test.csv'%(samplename),sep='\t',header=0)
    pd.DataFrame(PFCi.index).to_csv('%s_doublet_score_cell_name_test.csv'%(samplename),sep='\t',header=0)

