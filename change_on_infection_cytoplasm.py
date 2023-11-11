#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as ss

f=open("OVERREPRESENTATION_SPLASH/unique_overrepresented.csv").readlines()
overrepresented_genes=[]
for i in f:
    l=i.strip()
    overrepresented_genes.append(l)

f=open("Processing/gene_map.csv").readlines()
gene_map={}
for i in f:
    l=i.strip().split(",")
    gene_map[l[1]]=l[0]

#print(gene_map)

overrepresented_ids=[]
for i in overrepresented_genes:
    if i in gene_map:
        overrepresented_ids.append(gene_map[i])
    else:
        print("%s not found"%(i))

#print(len(overrepresented_genes))

f=open("Processing/res_contrast_infected-cytoplasm_vs_non_infected-cytoplasm.csv").readlines()
diffExp_cytoplasm={}
for i in f[1:]:
    l=i.strip().split()
    gname=l[0].replace("\"","")
    if float(l[1]) > 1000.0:
        diffExp_cytoplasm[gname]={'basemag':float(l[1]),'lfc':float(l[2]),'abs_lfc':abs(float(l[2]))}
        try:
            adj_pval=float(l[-1])
        except:
            adj_pval=1.0
        diffExp_cytoplasm[gname]['adj_pval']=-np.log10(adj_pval)

overrepresented_cytoplasm_basemag=[]
overrepresented_cytoplasm_abs_lfc=[]
overrepresented_cytoplasm_adj_pval=[]
nonrep_cytoplasm_basemag=[]
nonrep_cytoplasm_abs_lfc=[]
nonrep_cytoplasm_adj_pval=[]
for i in diffExp_cytoplasm:
    if i in overrepresented_ids:
        overrepresented_cytoplasm_basemag.append(diffExp_cytoplasm[i]['basemag'])
        overrepresented_cytoplasm_abs_lfc.append(diffExp_cytoplasm[i]['abs_lfc'])
        overrepresented_cytoplasm_adj_pval.append(diffExp_cytoplasm[i]['adj_pval'])
    else:
        nonrep_cytoplasm_basemag.append(diffExp_cytoplasm[i]['basemag'])
        nonrep_cytoplasm_abs_lfc.append(diffExp_cytoplasm[i]['abs_lfc'])
        nonrep_cytoplasm_adj_pval.append(diffExp_cytoplasm[i]['adj_pval'])

#print(overrepresented_cytoplasm_basemag)

t,p_mag=ss.ranksums(overrepresented_cytoplasm_basemag,nonrep_cytoplasm_basemag)
print("basemag: t: %f p: %f"%(t,p_mag))
t,p_lfc=ss.ranksums(overrepresented_cytoplasm_abs_lfc,nonrep_cytoplasm_abs_lfc)
print("abs_lfc: t: %f p: %f"%(t,p_lfc))
t,p_pval=ss.ranksums(overrepresented_cytoplasm_adj_pval,nonrep_cytoplasm_adj_pval)
print("adj_pval: t: %f p: %f"%(t,p_pval))

matplotlib.rcParams.update({
    'savefig.bbox':'tight',
    'font.size':16,
    'axes.linewidth':2,
    'xtick.major.width':2,
    'ytick.major.width':2,
})

plt.figure(1,figsize=(12,12))
plt.subplot2grid((1,3),(0,0))
plt.title("p=%g"%(p_mag))
plt.violinplot([
    overrepresented_cytoplasm_basemag,
    nonrep_cytoplasm_basemag,
            ],showmedians=True)
plt.xticks([1,2],[
    "Base Expression\n(interactors)",
    "Base Expression\n(non-interactors)",
],rotation='vertical')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.subplot2grid((1,3),(0,1))
plt.title("Cytoplasm\np=%g"%(p_lfc))
plt.violinplot([
    overrepresented_cytoplasm_abs_lfc,
    nonrep_cytoplasm_abs_lfc,
            ],showmedians=True)
plt.xticks([1,2],[
    "Absolute LogFC\n(interactors)",
    "Absolute LogFC\n(non-interactors)",
],rotation='vertical')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.subplot2grid((1,3),(0,2))
plt.title("p=%g"%(p_pval))
plt.violinplot([
    overrepresented_cytoplasm_adj_pval,
    nonrep_cytoplasm_adj_pval,
            ],showmedians=True)
plt.xticks([1,2],[
    "-log10(p-value)\n(interactors)",
    "-log10(p-value)\n(non-interactors)",
],rotation='vertical')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.savefig("Overrepresented_gene_expression_cytoplasm.pdf",format='pdf',dpi=300)

#print(diffExp_cytoplasm)
