import pandas as pd
import sys

dataFile = sys.argv[1]

# Import data file
data = pd.read_csv(dataFile, sep='\t')

data = data.iloc[:, 2:] 
column_names = data.columns.tolist()

output_file = sys.argv[2]

with open(output_file, 'w') as f:
    for column in column_names:
        f.write(column + '\n')