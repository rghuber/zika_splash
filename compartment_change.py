#!/usr/bin/env python3
import sys
import copy
import numpy as np
import scipy.stats as ss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    'savefig.bbox':'tight',
    'font.size':20,
})

f=open("ZikvGeneAbundance_TPM.csv").readlines()

compartment=[
    'cytoplasm',
    'cytoplasm',
    'nucleoplasm',
    'nucleoplasm',
    'chromatin',
    'chromatin',
    'full_lysate',
    'full_lysate',
    'cytoplasm',
    'cytoplasm',
    'nucleoplasm',
    'nucleoplasm',
    'chromatin',
    'chromatin',
    'full_lysate',
    'full_lysate',
]

infected=[
    'infected',
    'infected',
    'infected',
    'infected',
    'infected',
    'infected',
    'infected',
    'infected',
    'non_infected',
    'non_infected',
    'non_infected',
    'non_infected',
    'non_infected',
    'non_infected',
    'non_infected',
    'non_infected',
]

label_list=[
    'cytoplasm',
#    'cytoplasm 2',
    'nucleoplasm',
#    'nucleoplasm 2',
    'chromatin',
#    'chromatin 2',
]

expression_table={
    'name':"",
    'infected':
        {
        'cytoplasm':[],
        'nucleoplasm':[],
        'chromatin':[],
        'full_lysate':[],
        },
    'non_infected':
        {
        'cytoplasm':[],
        'nucleoplasm':[],
        'chromatin':[],
        'full_lysate':[],
        }
}

data_table={}

gene_ids=[]

of=open("compartment_change.csv","w")
of.write("gene_id,gene_name,rmsd,normed_rmsd")
of.write(",,infected_cytoplasm,infected_nucleoplasm,infected_chromatin")
of.write(",,non_infected_cytoplasm,non_infected_nucleoplasm,non_infected_chromatin")
#for n in range(len(compartment)):
#    if n%8 == 0:
#        of.write(",")
#    if compartment[n] != 'full_lysate':
#        of.write(",%s-%s"%(infected[n],compartment[n]))
of.write("\n")

def rmsd(a,b):
    rmsd=0.0
    for i in range(len(a)):
        rmsd+=(a[i]-b[i])**2.0
    return (rmsd**0.5)/float(len(a))

def normed_rmsd(a,b):
    rmsd=0.0
    a=np.array(a)/np.linalg.norm(a)
    b=np.array(b)/np.linalg.norm(b)
    for i in range(len(a)):
        rmsd+=(a[i]-b[i])**2.0
    return (rmsd**0.5)/float(len(a))


for i in f[1:]:
    l=i.split(",")
    gene_id=l[0]
    gene_ids.append(gene_id)
    gene_name=l[1]
    data_table[gene_id]=copy.deepcopy(expression_table)
    data_table[gene_id]['name']=gene_name
    values=[float(n) for n in l[2:]]
    
    for n in range(len(values)):
        data_table[gene_id][infected[n]][compartment[n]].append(values[n])
    for n in ['infected','non_infected']:
        for m in ['cytoplasm','nucleoplasm','chromatin','full_lysate']:
            data_table[gene_id][n][m]=np.mean(data_table[gene_id][n][m])
    
    t_infected=[
        data_table[gene_id]['infected']['cytoplasm'],
        data_table[gene_id]['infected']['nucleoplasm'],
        data_table[gene_id]['infected']['chromatin'],
    #    data_table[gene_id]['infected']['cytoplasm'][0],
    #    data_table[gene_id]['infected']['cytoplasm'][1],
    #    data_table[gene_id]['infected']['nucleoplasm'][0],
    #    data_table[gene_id]['infected']['nucleoplasm'][1],
    #    data_table[gene_id]['infected']['chromatin'][0],
    #    data_table[gene_id]['infected']['chromatin'][1],
    #    data_table[gene_id]['infected']['full_lysate'][0],
    #    data_table[gene_id]['infected']['full_lysate'][1],
    ]
    
    t_non_infected=[
        data_table[gene_id]['non_infected']['cytoplasm'],
        data_table[gene_id]['non_infected']['nucleoplasm'],
        data_table[gene_id]['non_infected']['chromatin'],
    #    data_table[gene_id]['non_infected']['cytoplasm'][0],
    #    data_table[gene_id]['non_infected']['cytoplasm'][1],
    #    data_table[gene_id]['non_infected']['nucleoplasm'][0],
    #    data_table[gene_id]['non_infected']['nucleoplasm'][1],
    #    data_table[gene_id]['non_infected']['chromatin'][0],
    #    data_table[gene_id]['non_infected']['chromatin'][1],
    #    data_table[gene_id]['non_infected']['full_lysate'][0],
    #    data_table[gene_id]['non_infected']['full_lysate'][1],
    ]

    #normalization
    #t_infected=np.array(t_infected)/np.linalg.norm(t_infected)
    #t_non_infected=np.array(t_non_infected)/np.linalg.norm(t_non_infected)

    if sum(t_infected) < 1.0 or sum(t_non_infected) < 1.0:
        continue

    #r,p=ss.pearsonr(t_infected,t_non_infected)
    rmsd_v=rmsd(t_infected,t_non_infected)
    normed_rmsd_v=normed_rmsd(t_infected,t_non_infected)

    t_infected=np.array(t_infected)
    t_non_infected=np.array(t_non_infected)

    if rmsd_v > 5.0 and normed_rmsd_v > 0.15 and not (any(t_infected==0) or any(t_non_infected==0)):
        of.write("%s,%s,%.3f,%.3f"%(gene_id,gene_name,rmsd_v,normed_rmsd_v))
        of.write(",")
        for n in range(3):
            of.write(",%f"%(t_infected[n]))
        of.write(",")
        for n in range(3):
            of.write(",%f"%(t_non_infected[n]))
        of.write("\n")

        plt.figure(1,figsize=(12,12))
        maxmax=max([max(t_infected),max(t_non_infected)])*1.1
        plt.scatter(t_non_infected,t_infected,s=25,c='k')
        plt.plot([0,maxmax],[0,maxmax],ls='--',color='grey')
        plt.title("%s (%s) rmsd %.3f normed_rmsd %.3f"%(gene_name,gene_id,rmsd_v,normed_rmsd_v))
        plt.xlabel("TPM non-infected")
        plt.ylabel("TPM infected")
        plt.xlim(0,maxmax)
        plt.ylim(0,maxmax)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        for n in (range(3)):
            plt.annotate(label_list[n],xy=(t_non_infected[n]+maxmax*0.01,t_infected[n]+maxmax*0.01),fontsize='x-small')
        plt.savefig("CompartmentShift/Shift_%s.pdf"%(gene_id),format='pdf',dpi=300)
        plt.clf()