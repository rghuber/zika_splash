#!/usr/bin/env python3
import sys

gene_names={}
f=open("gene_map.csv").readlines()
for i in f:
    l=i.strip().split(",")
    gene_names[l[0]]=l[1]

gene_lengths={}
f=open("../GRCh38_INDEX/GRCh38_lengths.tsv").readlines()
for i in f[1:]:
    l=i.split()
    gene_lengths[l[0]]=float(l[2])/1000.0

dsets={}
counts={}

fns=[]

for i in range(292,308):
    fn="RHX%d"%(i)
    fns.append(fn)
    dsets[fn]={}
    counts[fn]=0
    f=open("RHX%d.Gene.out.tab"%(i)).readlines()
    for j in f[4:]:
        l=j.split()
        gene=l[0]
        if gene in gene_lengths:
            dsets[fn][gene]=float(l[1])/gene_lengths[gene]
            counts[fn]+=dsets[fn][gene]

genes=sorted(list(dsets[fn].keys()))

of=open("ZikvGeneAbundance_TPM.csv","w")
of.write("EnsemblID,GeneName")
for i in fns:
    of.write(",%s"%(i))
of.write("\n")

for n in counts:
    counts[n]=counts[n]/1.0e6


no_lengths=0
for i in genes:
    if i in gene_lengths:
        if i in gene_names:
            of.write("%s,%s"%(i,gene_names[i]))
        else:
            of.write("%s,-"%(i))
        for fn in fns:
            of.write(",%f"%(dsets[fn][i]/counts[fn]))
        of.write("\n")
    else:
        print("no length for %s"%(i))
        no_lengths+=1

print(no_lengths)

#    of.write()

#print(len(common))


