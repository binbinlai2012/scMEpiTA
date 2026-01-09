#!/usr/bin/env python

import os,sys
from string import *
import simplesam

inbam=sys.argv[1]
outsam=sys.argv[2]

barcode_tag='CB'
umi_tag='UB'

with simplesam.Reader(open(inbam, 'r')) as in_bam:
    with simplesam.Writer(open(outsam, 'w'), in_bam.header) as out_sam:
        for read in in_bam:
            read[umi_tag] = read.qname.split("_")[2] # add the umi tag
            read[barcode_tag] = read.qname.split("_")[1] # add the barcode tag
            out_sam.write(read)

