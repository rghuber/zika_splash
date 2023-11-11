#!/usr/bin/env python3
import sys

dsets=["RHX%03d"%(d) for d in range(292,308)]
print(dsets)


data={}

for ds in dsets:
    data[ds]={}
    f=open("%s.Gene.out.tab"%(ds)).readlines()
    for i in f:
        if "ENSG" in i:
            l=i.split()
            gene=l[0]
            count=int(l[3])
            data[ds][gene]=count

gene_table=sorted(list(data[ds].keys()))

of=open("gene_count_table.csv","w")
of.write("gene")
for ds in dsets:
    of.write(",%s"%(ds))
of.write("\n")
for gene in data[dsets[0]]:
    of.write("%s"%(gene))
    for ds in dsets:
        of.write(",%d"%(data[ds][gene]))
    of.write("\n")

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

ordered_data={}
for n,state in enumerate(infected):
    ordered_data[state]={n:[] for n in set(compartment)}
for n,state in enumerate(infected):
    ordered_data[state][compartment[n]].append(dsets[n])

        
contrasts=[
    [['infected','cytoplasm'],['non_infected','cytoplasm']],
    [['infected','nucleoplasm'],['non_infected','nucleoplasm']],
    [['infected','chromatin'],['non_infected','chromatin']],
    [['infected','full_lysate'],['non_infected','full_lysate']],

    [['infected','cytoplasm'],['infected','nucleoplasm']],
    [['infected','cytoplasm'],['infected','chromatin']],
    [['infected','cytoplasm'],['infected','full_lysate']],
    [['infected','nucleoplasm'],['infected','full_lysate']],
    [['infected','nucleoplasm'],['infected','chromatin']],
    [['infected','chromatin'],['infected','full_lysate']],
    
    [['non_infected','cytoplasm'],['non_infected','full_lysate']],
    [['non_infected','cytoplasm'],['non_infected','chromatin']],
    [['non_infected','cytoplasm'],['non_infected','full_lysate']],
    [['non_infected','nucleoplasm'],['non_infected','full_lysate']],
    [['non_infected','nucleoplasm'],['non_infected','chromatin']],
    [['non_infected','chromatin'],['non_infected','full_lysate']],
]

for ctrst in contrasts:
    fname="%s-%s_vs_%s-%s"%(ctrst[0][0],ctrst[0][1],ctrst[1][0],ctrst[1][1])
    of=open("gene_counts_%s.csv"%(fname),"w")
    of.write("gene")
    
    dset_list=[]
    for n in ctrst:
        dset_list+=ordered_data[n[0]][n[1]]
    for n in dset_list:
        of.write(",%s"%(n))
    of.write("\n")
    for gene in gene_table:
        of.write("%s"%(gene))
        for n in dset_list:
            of.write(",%d"%(data[n][gene]))
        of.write("\n")
    of.close()
    of=open("annotations_%s.csv"%(fname),"w")
    of.write("dset,contrast\n")
    of.write("%s,A\n"%(dset_list[0]))
    of.write("%s,A\n"%(dset_list[1]))
    of.write("%s,B\n"%(dset_list[2]))
    of.write("%s,B\n"%(dset_list[3]))
    of.close()



of=open("annotations.csv","w")
of.write("dset,compartment,zika_infected\n")
for i in range(len(dsets)):
    of.write("%s,%s,%s\n"%(dsets[i],compartment[i],infected[i]))