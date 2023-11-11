#!/usr/bin/env python3
import sys

f=open("gene_map.csv").readlines()

gene_type={}
for i in f:
    l=i.strip().split(",")
    gene_type[l[1]]=l[2]

f=open(sys.argv[1]).readlines()
for i in f:
    gname=i.strip()
    if gname in gene_type:
        print(gname,gene_type[gname])
    else:
        print(gname,"UNKNOWN")
