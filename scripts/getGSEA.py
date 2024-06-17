

import pandas as pd
import gseapy as gp
import yaml
import sys

# Sample list argument
if len(sys.argv) < 2 :
    print("Defaulting to all samples")
else :
    sampleList = sys.argv[1]

# Set up yaml
with open('./input/config.yaml', 'r') as f:
    config = yaml.safe_load(f)

# Import data file
file = pd.read_csv(config['dataFile'], sep='\t')
cleanfile = file.dropna()

# Import gene set file
with open(config['geneSetFile'], 'r') as file:
    geneSets = [line.strip() for line in file]

# Get list of samples
if len(sys.argv) == 2   :
    with open(sampleList, 'r') as file:
        column_names = [line.strip() for line in file]
else :
    column_names = cleanfile.columns.tolist()[2:]

cleanfile['Average'] = cleanfile[column_names].mean(axis=1)


# Get the dataframe for one sample
cleanfile_sample = {'Gene':cleanfile['Hugo_Symbol'], 'Sample':cleanfile["Average"] }
cleanfile_sample = pd.DataFrame(cleanfile_sample)


cleanfile_sample_sorted = cleanfile_sample.sort_values(by='Sample', ascending=False)

# Rank
pre_res = gp.prerank(rnk=cleanfile_sample_sorted, # or rnk = rnk,
                     gene_sets=geneSets,
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

outlist = outlist.sort_values(by='pval')

outlist.to_csv(config['resultFile'], index=False)
