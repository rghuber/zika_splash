#!/usr/bin/env python3

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

dsets={
    'DEN1':'EU081230',
    'ZAFR':'AY632535',
    'ZBRA':'KU497555',
    'ZILM':'KJ776791',
    'ZSIN':'KY241786',
}



for ds in dsets:
    f=open("%s.chimera-flt.csv"%(ds)).readlines()
    loc=np.zeros(10000)
    for ln in f:
        l=ln.strip().split()
        gene_a=l[2]
        start_a=int(l[3])
        stop_a=int(l[4])
        gene_b=l[9]
        start_b=int(l[10])
        stop_b=int(l[11])

        if gene_a == dsets[ds] and gene_b == "MALAT1":
            #print(start_b,stop_b)
            for n in range(start_b,stop_b):
                loc[n]+=1
        if gene_b == dsets[ds] and gene_a == "MALAT1":
            #print(start_a,stop_a,l)
            for n in range(start_a,stop_a):
                loc[n]+=1

    print(ds,max(loc[8200:8300]))

    plt.figure(1,figsize=(20,4))
    plt.title("%s-MALAT1 on MALAT1"%(ds))
    plt.subplot2grid((3,1),(0,0),rowspan=2)
    plt.plot(range(len(loc)),loc,color='k')
    plt.xlim(0,9000)
    plt.xticks(range(0,9001,1000),[])
    plt.ylabel("SPLASH")
    plt.subplot2grid((3,1),(2,0))
    plt.bar([8200],[1],width=100,color='r')
    plt.xticks(range(0,9001,1000))
    plt.yticks([0,1],[])
    plt.xlim(0,9000)
    plt.xlabel("MALAT1 coordinate")
    plt.savefig("MALAT1_location_%s.pdf"%(ds),format='pdf',dpi=300)
    plt.clf()

