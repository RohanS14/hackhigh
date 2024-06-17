############################################################
#                                                          #
#   Hackhigh - GSEA analysis from cbioportal mRNA samples  #
#                                                          #
############################################################

How to run on end to end script:

1. Edit config.yaml (/scripts/input/config.yaml)

   - Path to cbioportal data
   - Path to list of gene sets to use
   - Path to output file

2. Run!

    - 'cd scripts'
    - 'python getGSEA.py {path to sample list}'
    - By default, the program will use all samples in the data file. Provide a list of samples to select samples within the data file.

How to run gsea.ipynb notebook:

1. Edit dataFile, geneSetFile, resultFile

2. Run 'Run GSEA'

3. Run 'plotGSEA({Gene Set Term})' for GSEA enrichment plot for one gene set

How to run visualize.ipynb notebook:

1. First run gsea.ipynb or the end to end script

2. Run the script (:
