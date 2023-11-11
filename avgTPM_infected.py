#!/usr/bin/env python3
import sys
import numpy as np

f=open("ZikvGeneAbundance_TPM.csv").readlines()
of=open("InfectedGeneAbundance.csv","w")
of.write("gene_id,gene_name,cytoplasm,nucleoplasm,chromatin,full_lysate\n")
for i in f[1:]:
    l=i.strip().split(",")
    of.write("%s,%s"%(l[0],l[1]))
    of.write(",%f"%(np.mean([float(l[2]),float(l[3])])))
    of.write(",%f"%(np.mean([float(l[4]),float(l[5])])))
    of.write(",%f"%(np.mean([float(l[6]),float(l[7])])))
    of.write(",%f"%(np.mean([float(l[8]),float(l[9])])))
    of.write("\n")
