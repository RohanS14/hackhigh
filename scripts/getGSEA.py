

import pandas as pd
import gseapy as gp


# Import file
file = pd.read_csv("../data/data_mrna_seq_v2_rsem_zscores_ref_diploid_samples.txt", sep='\t')
cleanfile = file.dropna()

# Get the dataframe for one sample
cleanfile_sample = {'Gene':cleanfile['Hugo_Symbol'], 'Sample':cleanfile["TCGA-OR-A5J1-01"] }
cleanfile_sample = pd.DataFrame(cleanfile_sample)

cleanfile_sample_sorted = cleanfile_sample.sort_values(by='Sample', ascending=False)

# Rank
pre_res = gp.prerank(rnk=cleanfile_sample_sorted, # or rnk = rnk,
                     gene_sets='GO_Molecular_Function_2023',
                     seed=6
                    )

# Results
outlist=[]
for term in pre_res.results.keys():
    p=pre_res.results[term]['pval']
    fdr=pre_res.results[term]['fdr']
    nes=pre_res.results[term]['nes']
    es=pre_res.results[term]['es']
    gene=pre_res.results[term]['lead_genes']
    outlist.append([term,p,fdr,nes,es,gene])

outlist=pd.DataFrame(outlist, columns=['Term', 'pval', 'fdr','nes','es','gene'])

outlist.sort_values(by='pval')

outlist.to_csv('../results/gsea1.csv', index=False)