#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import venn



f=open("InfectedGeneAbundance.csv").readlines()
compartments=f[0].strip().split(",")[2:]
gene_abundance={n:{} for n in compartments}
gene_abundance['gene_id']={}

gene_names=[]

for i in f[1:10]:
    l=i.strip().split(",")
    gene_id=l[0]
    gene_name=l[1]
    gene_names.append(gene_name)
    gene_rpkm=[float(n) for n in l[2:]]
    gene_abundance['gene_id'][gene_name]=gene_id
    for n,c in enumerate(compartments):
        gene_abundance[c][gene_name]=gene_rpkm[n]

viruses=[
    'DEN1',
    'ZAFR',
    'ZBRA',
    'ZILM',
    'ZSIN',
]

overrepresented_genes={n:{m:[] for m in viruses} for n in compartments}
overrepresented_counts={n:{m:[] for m in viruses} for n in compartments}
overrepresented_peaks={n:{m:[] for m in viruses} for n in compartments}

for virus in viruses:
    for compartment in compartments:
        f=open("overrepresented_genes_%s-%s.csv"%(virus,compartment)).readlines()
        for i in f[1:]:
            l=i.strip().split(",")
            overrepresented_genes[compartment][virus].append(l[0])
            overrepresented_counts[compartment][virus].append(int(l[1]))
            overrepresented_peaks[compartment][virus].append(int(l[3]))

matplotlib.rcParams.update({
    'savefig.bbox':'tight',
})


mesh_plot_tables={n:[] for n in compartments}
for nc,compartment in enumerate(compartments):
    total_set=set(overrepresented_genes[compartment][viruses[0]])
    for n in range(1,len(viruses)):
        total_set=total_set.union(set(overrepresented_genes[compartment][viruses[n]]))
    total_set=list(total_set)
    for gene in total_set:
        mesh_plot_tables[compartment].append([])
        for n in range(len(viruses)):
            if gene in overrepresented_genes[compartment][viruses[n]]:
                mesh_plot_tables[compartment][-1].append(1)
            else:
                mesh_plot_tables[compartment][-1].append(0)
    sortorder_input=[]
    for i in range(len(total_set)):
        sum5=sum(mesh_plot_tables[compartment][i])
        sum4=sum(mesh_plot_tables[compartment][i][1:])
        sum3=sum(mesh_plot_tables[compartment][i][:3])
        sum2=sum(mesh_plot_tables[compartment][i][:2])
        sortorder_input.append(1000*sum5+100*sum4+10*sum3+sum2)
    sortorder_input=np.argsort(sortorder_input)
    sorted_set=[total_set[n] for n in sortorder_input]
    sorted_table=[mesh_plot_tables[compartment][n] for n in sortorder_input]
    mesh_plot_tables[compartment]=sorted_table

    plt.figure(1,figsize=(6,48))

    plt.pcolormesh(mesh_plot_tables[compartment],cmap=plt.cm.Reds,vmin=0,vmax=1)
    plt.xticks(np.linspace(0.5,len(viruses)-0.5,len(viruses)),viruses)
    plt.yticks(np.linspace(0.5,len(sorted_set)-0.5,len(sorted_set)),sorted_set,fontsize='x-small')

    for i in range(len(viruses)):
        plt.plot([i,i],[0,len(sorted_set)],color='k')
    for i in range(len(sorted_set)):
        plt.plot([0,len(viruses)],[i,i],color='k')

    plt.savefig("mesh_plot_%s.pdf"%(compartment),format='pdf',dpi=300)
    plt.savefig("mesh_plot_%s.svg"%(compartment),format='pdf',dpi=300)

    plt.clf()
