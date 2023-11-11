#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

gene_map={}
f=open("gene_map.csv").readlines()
for i in f:
    l=i.strip().split(",")
    gene_map[l[0]]=l[1]

data={}

f=open(sys.argv[1]).readlines()
of=open("ordered_%s"%(sys.argv[1]),"w")

x=[]
y=[]

of.write("gene_id,gene_name,basemag,lfc,p,padj\n")
for i in f[1:]:
    l=i.strip().split()
    gene=l[0].replace("\"","")
    if gene in gene_map:
        name=gene_map[gene]
    else:
        name=gene
    basemag=float(l[1])
    lfc=float(l[2])
    try:
        p=float(l[-2])
    except:
        p=1.0
    try:
        padj=float(l[-1])
    except:
        padj=1.0
    if abs(lfc) > 1.0 and padj < 0.05:
        of.write("%s,%s,%f,%f,%g,%g\n"%(gene,name,basemag,lfc,p,padj))
    x.append(lfc)
    if padj > 0:
        log_p=-np.log10(padj)
    else:
        log_p=250
    y.append(log_p)

y_scale=max(y)*1.1
x_scale=max([abs(max(x)),abs(min(x))])*1.1

plt.figure(1,figsize=(12,12))
plt.scatter(x,y,s=5,c='k')
plt.plot([-1,-1],[0,y_scale],ls='--',color='grey')
plt.plot([1,1],[0,y_scale],ls='--',color='grey')
plt.plot([-x_scale,x_scale],[-np.log10(0.05),-np.log10(0.05)],ls='--',color='grey')
plt.xlim(-x_scale,x_scale)
plt.ylim(0,y_scale)
plt.savefig("volcano_%s.pdf"%(sys.argv[1].replace(".csv","")),format='pdf',dpi=300)