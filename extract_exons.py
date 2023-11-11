#!/usr/bin/env python3
import sys
import subprocess as sp

f=open("exons.tsv").readlines()
for i in f:
    l=i.strip().split()
    gene=l[0]
    chr=l[1].split("_")[0]
    strand=l[2]
    starts=[int(n) for n in l[3].strip(",").split(",")]
    stops=[int(n) for n in l[4].strip(",").split(",")]
    exon_cons=[]
    ex_len=0
    samp_len=0
    for i in range(len(starts)):
        sp.call("bigWigToBedGraph -chrom=%s -start=%d -end=%d ~/ZikvPhastCons/hg19.100way.phastCons.bw temp.bed"%(chr,starts[i],stops[i]),shell=True)
        f=open("temp.bed").readlines()
        vals=[]
        for j in f:
            l=j.split()
            bedstart=int(l[1])
            bedstop=int(l[2])
            val=float(l[3])
            for k in range(bedstart,bedstop):
                vals.append(val)
        exon_cons.append(vals)
        #print(chr,len(vals),stops[i]-starts[i])
        if len(vals) != stops[i] - starts[i]:
            print("issue with %s %s %s %s"%(gene,chr,starts[i],stops[i]))
        ex_len+=len(vals)
        samp_len+=stops[i]-starts[i]

    final_cons=[]
    if strand == "+":
        for n in exon_cons:
            final_cons=final_cons+n
    else:
        for n in reversed(exon_cons):
            final_cons=final_cons+list(reversed(n))
    if samp_len==ex_len:
        of=open("phastCons_%s.csv"%(gene),"w")
        for n,val in enumerate(final_cons):
            of.write("%s,%f\n"%(n+1,val))
    print(gene,ex_len,samp_len,len(final_cons))


