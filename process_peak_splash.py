#!/usr/bin/env python3 
import sys
import numpy as np

dsets={
    "DEN1":"EU081230",
    "ZAFR":"AY632535",
    "ZBRA":"KU497555",
    "ZILM":"KJ776791",
    "ZSIN":"KY241786",
}

dlocs={ds:{} for ds in dsets}
glocs={ds:{} for ds in dsets}
gcnt={ds:{} for ds in dsets}

dmax={ds:{} for ds in dsets}
gmax={ds:{} for ds in dsets}

f=open("gene_map.csv").readlines()
gene_type={}
for i in f:
    l=i.strip().split(",")
    gene_type[l[1]]=l[2]

for ds in dsets:
    f=open("%s.chimera-flt.csv"%(ds)).readlines()
    for i in f:
        l=i.strip().split()
        gene_a=l[2]
        start_a=int(l[3])
        stop_a=int(l[4])
        gene_b=l[9]
        start_b=int(l[10])
        stop_b=int(l[11])
        if (gene_a == dsets[ds]) and (gene_b != dsets[ds]):
            if gene_b not in dlocs[ds]:
                gcnt[ds][gene_b]=1
                dlocs[ds][gene_b]=[0 for i in range(11000)]
                glocs[ds][gene_b]=[0 for i in range(100000)]
                for n in range(start_a,stop_a):
                    dlocs[ds][gene_b][n]+=1
                for n in range(start_b,stop_b):
                    glocs[ds][gene_b][n]+=1
            else:
                gcnt[ds][gene_b]+=1
                for n in range(start_a,stop_a):
                    dlocs[ds][gene_b][n]+=1
                for n in range(start_b,stop_b):
                    glocs[ds][gene_b][n]+=1

        if (gene_b == dsets[ds]) and (gene_a != dsets[ds]):
            if gene_a not in dlocs[ds]:
                gcnt[ds][gene_a]=1
                dlocs[ds][gene_a]=[0 for i in range(11000)]
                glocs[ds][gene_a]=[0 for i in range(100000)]
                for n in range(start_b,stop_b):
                    dlocs[ds][gene_a][n]+=1
                for n in range(start_a,stop_a):
                    glocs[ds][gene_a][n]+=1
            else:
                gcnt[ds][gene_a]+=1
                for n in range(start_b,stop_b):
                    dlocs[ds][gene_a][n]+=1
                for n in range(start_a,stop_a):
                    glocs[ds][gene_a][n]+=1

    for gene in dlocs[ds]:
        maxloc=np.argmax(dlocs[ds][gene])
        maxval=dlocs[ds][gene][maxloc]
        dmax[ds][gene]={'highest_peak':maxval,'peak_location':maxloc}
    for gene in glocs[ds]:
        maxloc=np.argmax(glocs[ds][gene])
        maxval=glocs[ds][gene][maxloc]
        gmax[ds][gene]={'highest_peak':maxval,'peak_location':maxloc}
    
    of=open("peak_locations_%s.csv"%(ds),"w")
    of.write("gene,total_count,denzik_peak_height,denzik_peak_location,host_peak_height,host_peak_location,gene_type\n")
    for gene in dmax[ds]:
        of.write("%s,%d,%d,%d,%d,%d"%(gene,gcnt[ds][gene],dmax[ds][gene]['highest_peak'],dmax[ds][gene]['peak_location'],gmax[ds][gene]['highest_peak'],gmax[ds][gene]['peak_location']))
        if gene in gene_type:
            of.write(",%s\n"%(gene_type[gene]))
        else:
            of.write(",\n")