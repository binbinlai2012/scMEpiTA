#!/usr/bin/env python

import os,sys
import loompy
import numpy as np

loomfile=sys.argv[1]
outfile=sys.argv[2]

ds = loompy.connect(loomfile)
(umic,) = ds.map([np.sum], axis=1)
cell_ids = ds.ca['CellID']
cell_umi_dict = dict(zip(cell_ids, umic))
sorted_cell_umi = dict(sorted(cell_umi_dict.items(), key=lambda x: x[1], reverse=True))

(splicedc,)=ds.layers["spliced"].map([np.sum], axis=1)
(unsplicedc, ) = ds.layers['unspliced'].map([np.sum], axis=1)
(ambiguousc, ) = ds.layers['ambiguous'].map([np.sum], axis=1)

cell_splicedc_dict = dict(zip(cell_ids, splicedc))
cell_unsplicedc_dict = dict(zip(cell_ids, unsplicedc))
cell_ambiguousc_dict = dict(zip(cell_ids, ambiguousc))

outf = open(outfile, 'w')
outf.write("CellID\tUMI\tSpiced\tUnspliced\tAmbiguous\n")
for i, (cell_id, umi_count) in enumerate(list(sorted_cell_umi.items())):
    outf.write("%s\t%d\t%d\t%d\t%d\n"%(cell_id, umi_count, 
                                           cell_splicedc_dict[cell_id], 
                                           cell_unsplicedc_dict[cell_id], 
                                           cell_ambiguousc_dict[cell_id]))
outf.close()
ds.close()
