#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import venn

viruses=[
    'DEN1',
    'ZAFR',
    'ZBRA',
    'ZILM',
    'ZSIN',
]

compartments=[
    'cytoplasm',
    'nucleoplasm',
    'chromatin',
    'full_lysate',
]

gene_lists={}
for virus in viruses:
    gene_lists[virus]={}
    for compartment in compartments:
        gene_lists[virus][compartment]=[]
        f=open("overrepresented_splash_%s_%s.csv"%(virus,compartment)).readlines()
        for i in f[1:]:
            l=i.split(",")
            gene_lists[virus][compartment].append(l[0])
        gene_lists[virus][compartment]=set(gene_lists[virus][compartment])
        print(virus,compartment,len(gene_lists[virus][compartment]))

genes_by_compartment={n:{m:gene_lists[m][n] for m in viruses} for n in compartments}


plt.figure(1,figsize=(36,12))
all_venns=[]
for n,compartment in enumerate(compartments):
    plt.subplot2grid((2,5),(0,n))
    plt.title("overrepresentation\n%s"%(compartment))
    ax=plt.gca()
    venn.venn(genes_by_compartment[compartment],ax=ax)

    center_venn=genes_by_compartment[compartment][viruses[0]]
    for i in range(1,len(viruses)):
        center_venn=center_venn.intersection(genes_by_compartment[compartment][viruses[i]])

    of=open("center_venn_peak_%s.txt"%(compartment),"w")
    for i in sorted(center_venn):
        of.write("%s\n"%(i))
    all_venns.append(center_venn)

all_center_venns=set(all_venns[0])
for i in range(1,len(all_venns)):
    print(len(all_venns[i]))
    all_center_venns=all_center_venns.intersection(all_venns[i])
of=open("all_center_venn_overlap_compartments_peak.txt","w")
for i in sorted(all_center_venns):
        of.write("%s\n"%(i))

all_venns=[]
for n,virus in enumerate(viruses):
    plt.subplot2grid((2,5),(1,n))
    plt.title("overrepresentation\n%s"%(virus))
    ax=plt.gca()
    venn.venn(gene_lists[virus],ax=ax)
    center_venn=gene_lists[virus][compartments[0]]
    for i in range(1,len(compartments)):
        center_venn=center_venn.intersection(gene_lists[virus][compartments[i]])
    of=open("center_venn_peak_%s.txt"%(virus),"w")
    for i in sorted(center_venn):
        of.write("%s\n"%(i))
    all_venns.append(center_venn)

all_center_venns=set(all_venns[0])
for i in range(1,len(all_venns)):
    print(len(all_venns[i]))
    all_center_venns=all_center_venns.intersection(all_venns[i])
of=open("all_center_venn_overlap_viruses_peak.txt","w")
for i in sorted(all_center_venns):
        of.write("%s\n"%(i))

plt.savefig("venn_overrepresentation_peak.pdf",format='pdf',dpi=300)