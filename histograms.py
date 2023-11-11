#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

t_vals=[float(n) for n in open("t_stat.csv").readlines()]
p_vals=[float(n) for n in open("p_stat.csv").readlines()]

neg_log_p=[]
for val in p_vals:
    if val != 0:
        neg_log_p.append(-np.log10(val))

y,x=np.histogram(t_vals,range=(-15,15),bins=30)
xlabels=[]
for i in range(len(x)-1):
    xlabels.append("%.2f-%.2f"%(x[i],x[i+1]))
plt.figure(1,figsize=(12,12))
plt.title("t values")
plt.bar(range(len(y)),y,color='k',width=0.8)
plt.xticks(range(len(y)),xlabels,rotation='vertical')
plt.ylabel("Frequency")
plt.savefig("t_value_distribution.pdf",format='pdf',dpi=300)

y,x=np.histogram(neg_log_p,range=(0,150),bins=30)
xlabels=[]
for i in range(len(x)-1):
    xlabels.append("%.2f-%.2f"%(x[i],x[i+1]))
plt.figure(1,figsize=(12,12))
plt.title("-log10(p) values")
plt.bar(range(len(y)),y,color='k',width=0.8)
plt.xticks(range(len(y)),xlabels,rotation='vertical')
plt.ylabel("Frequency")
plt.savefig("p_value_distribution.pdf",format='pdf',dpi=300)
