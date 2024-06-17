import pandas as pd
import gseapy as gp
import yaml
import sys

def getGSEA(geneSetFile, dataFile, sampleList=None, resultFile=None):
    
    # Get gene sets
    with open(geneSetFile, 'r') as file:
        geneSets = [line.strip() for line in file]
        
    # Get list of samples
    dataFile = pd.read_csv(dataFile, sep='\t')
    
    cleanfile = dataFile.dropna()

    if sampleList is not None:
        with open(sampleList, 'r') as file:
            column_names = [line.strip() for line in file]
    else:
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

    if resultFile is not None:
        outlist.to_csv(resultFile, index=False)
    
    return pre_res

def main():
    # Sample list argument
    if len(sys.argv) < 2 :
        print("Defaulting to all samples")
    else :
        sampleList = sys.argv[1]

    # Set up yaml
    with open('./input/config.yaml', 'r') as f:
        config = yaml.safe_load(f)

    # Import data file
    dataFile = pd.read_csv(config['dataFile'], sep='\t')
    
    # Import gene set file
    geneSetFile = config['geneSetFile']
    
    # Result file
    result_file = config['resultFile']
    
    # Run GSEA
    pre_res = getGSEA(geneSetFile, dataFile, sampleList=None, resultFile=None)
    
if __name__ == "__main__":
    main()