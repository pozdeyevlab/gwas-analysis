import polars as pl 
from pathlib import Path

# just a simple helper script to make the map file based on the available files in a given directory
dfs = []

files = Path('.').glob('*biobanks*')
for file in files:
	dct = {'phenotype':f"{file.name.split('_genes_')[1].replace('.tsv', '')}", 'ancestry':f"{file.name.split('_biobanks_')[1].split('_')[0]}", 'sex':f"{file.name.split('_leave')[0].split('_')[-1]}", 'file':f"{file.resolve()}"}
	dfs.append(pl.DataFrame(dct))
map = pl.concat(dfs, how='vertical')
print(map)
map.write_csv('map.tsv', separator='\t')
