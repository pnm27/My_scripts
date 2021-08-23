#!/usr/bin/env python

import scanpy as sc, anndata as ad, os

# To produce list of files, one can use: os.listdir, os.path.isfile and os.path.join 
samples = ['sample1.h5ad', 'sample2.h5ad', 'sample3.h5ad', 'sample4.h5ad']

adatas = [ad.read(samples) for sample in samples] # create anndata array

adata = ad.concat(adatas[:], index_unique='-', label='batch', keys=sample_names, join='inner') # it will add a 'batch' variable with keys as its value and will inner join on the 'genes'

# If this line produces error (unable to write file) then make sure adata.var_names is not same as one of the vars column
adata.write('combined.h5ad')