#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import venn

dsets={
    'infected':{},
    'compartment':{},
}

dsets['compartment']=['cytoplasm_vs_whole','nucleoplasm_vs_whole','chromatin_vs_whole']
dsets['infected']=['cytoplasm','nucleoplasm','whole']

gene_map={}
f=open("gene_map.csv").readlines()
for i in f:
    l=i.strip().split(",")
    gene_map[l[0]]=l[1]

data={}
for ds in dsets:
    data[ds]={}
    for crit in dsets[ds]:
        f=open("res_%s_%s.csv"%(ds,crit)).readlines()
        data[ds][crit]={'gene':[],'lfc':[],'pval':[],'basemag':[]}
        for i in f[1:]:
            l=i.split()
            gene=l[0].replace("\"","")
            lfc=float(l[2])
            basemag=float(l[1])
            try:
                if float(l[-1]) > 0:
                    pval=-np.log10(float(l[-1]))
                else:
                    pval=100
            except:
                pval=0.0
            data[ds][crit]['gene'].append(gene)
            data[ds][crit]['lfc'].append(lfc)
            data[ds][crit]['pval'].append(pval)
            data[ds][crit]['basemag'].append(basemag)


for ds in dsets:
    for crit in dsets[ds]:
        g=data[ds][crit]['gene']
        x=data[ds][crit]['lfc']
        y=data[ds][crit]['pval']
        bm=data[ds][crit]['basemag']
        
        of=open("volcano_table_%s_%s.csv"%(ds,crit),"w")
        of.write("gene_id,gene_name,baseMag,logFC,-log10(pval)\n")
        for i in range(len(g)):
            if abs(x[i]) > 1 and y[i] > 1.3:
                if g[i] in gene_map:
                    of.write("%s,%s,%f,%f,%f\n"%(g[i],gene_map[g[i]],bm[i],x[i],y[i]))
                else:
                    of.write("%s,%s,%f,%f,%f\n"%(g[i],g[i],bm[i],x[i],y[i]))
        of.close()
        
        x_scale=max([abs(min(x)),abs(max(x))])*1.1
        y_scale=max(y)*1.1
        plt.figure(1,figsize=(12,12))
        plt.title("%s %s"%(ds,crit))
        plt.scatter(x,y,s=5,c='k')
        plt.plot([-x_scale,x_scale],[1.3,1.3],ls='--',color='grey')
        plt.plot([-1,-1],[0,y_scale],ls='--',color='grey')
        plt.plot([1,1],[0,y_scale],ls='--',color='grey')
        plt.xlabel(r"log$_2$-fold change")
        plt.ylabel(r"-log$_{10}$ p")
        plt.xlim(-x_scale,x_scale)
        plt.ylim(0,y_scale)
        plt.savefig("volcano_%s_%s.pdf"%(ds,crit),format='pdf',dpi=300)
        plt.clf()


deg={}
for ds in dsets:
    deg[ds]={}
    for crit in dsets[ds]:
        deg[ds][crit]={'gene':[],'lfc':[],'pval':[],'basemag':[]}
        for n in range(len(data[ds][crit]['gene'])):
                if data[ds][crit]['pval'][n] > 1.3 and abs(data[ds][crit]['lfc'][n]) > 1.0:
                    deg[ds][crit]['gene'].append(data[ds][crit]['gene'][n])
                    deg[ds][crit]['lfc'].append(data[ds][crit]['lfc'][n])
                    deg[ds][crit]['pval'].append(data[ds][crit]['pval'][n])
                    deg[ds][crit]['basemag'].append(data[ds][crit]['basemag'][n])


matplotlib.rcParams.update({
    'savefig.bbox':'tight',
    'font.size':20,  
})

plt.figure(2,figsize=(12,12))
plt.title("%s"%(ds))
ax=plt.gca()
deg_venn={}
for crit in deg['compartment']:
    deg_venn[crit]=set(deg['compartment'][crit]['gene'])
venn.venn(deg_venn,ax=ax)
plt.savefig("venn_by_compartment.pdf",format='pdf',dpi=300)
plt.clf()

plt.figure(3,figsize=(12,12))
plt.title("%s"%(ds))
ax=plt.gca()
deg_venn={}
for crit in deg['infected']:
    deg_venn[crit]=set(deg['infected'][crit]['gene'])
venn.venn(deg_venn,ax=ax)
plt.savefig("venn_by_infected.pdf",format='pdf',dpi=300)
plt.clf()
