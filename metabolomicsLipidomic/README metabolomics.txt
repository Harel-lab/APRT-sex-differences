raw data:
metabolomics: inputFiles/Polar compounds_norm_weight.csv
lipidomic: inputFiles/Lipids_norm_weight.csv 
sample information: inputFiles/sample factors.csv

The raw data of either metabolomics and lipidomic are normalized with internal standard and sample weight.
The metabolites and lipids with less expression than 70% (19/27) were filtered out. After, they are transforme by Log2 and normalized by the average of wach metabolite or lipids.

