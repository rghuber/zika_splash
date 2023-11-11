#!/usr/bin/env python3
import sys
import os
import numpy as np
import scipy.stats as ss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import scipy.ndimage as snd

matplotlib.rcParams.update({
    'savefig.bbox':'tight',
    'font.size':16,
    'axes.linewidth':2,
    'xtick.major.width':2,
    'ytick.major.width':2,
})

f=open("../unique_overrepresented.csv").readlines()
gene_list=[]
for i in f:
    gene_list.append(i.strip())

window=20

dsets={
    "DEN1":"EU081230",
    "ZAFR":"AY632535",
    "ZBRA":"KU497555",
    "ZILM":"KJ776791",
    "ZSIN":"KY241786",
}

f=open("hg19_genes.tsv").readlines()
gene_loc={}
for i in f[1:]:
    l=i.strip().split()
    gname=l[12]
    glen=int(l[4])
    txstart=int(l[4])
    txstop=int(l[5])
    gchr=l[2]
    exon_starts=l[9]
    exon_stops=l[10]
    gstrand=l[3]
    gene_loc[gname]={'chr':gchr,
                     'txstart':txstart,
                     'txstop':txstop,
                     'exon_starts':exon_starts,
                     'exon_stops':exon_stops,
                     'strand':gstrand,
                     'length':glen,
                     }


ds_keys=list(dsets.keys())

gene_host_locations={n:{k:[0 for n in range(11000)] for k in dsets} for n in gene_list}
gene_host_positions={n:{k:[0 for n in range(25000)] for k in dsets} for n in gene_list}

for ds in dsets:
    vname=dsets[ds]
    f=open("../%s.chimera-flt.csv"%(ds)).readlines()
    for i in f:
        l=i.split()
        gene_a=l[2]
        start_a=int(l[3])-1
        stop_a=int(l[4])-1
        gene_b=l[9]
        start_b=int(l[10])-1
        stop_b=int(l[11])-1
        if gene_a in gene_list and gene_b == vname:
            for n in range(start_b,stop_b):
                gene_host_locations[gene_a][ds][n]+=1
            for n in range(start_a,stop_a):
                gene_host_positions[gene_a][ds][n]+=1
        if gene_b in gene_list and gene_a == vname:
            for n in range(start_a,stop_a):
                gene_host_locations[gene_b][ds][n]+=1
            for n in range(start_b,stop_b):
                gene_host_positions[gene_b][ds][n]+=1

def smooth(seq):
    smooth_seq=[]
    for i in range(len(seq)-window):
        smooth_seq.append(np.mean(seq[i:i+window]))
    return np.array([0 for n in range(window//2)]+smooth_seq+[0 for n in range(window//2)])


def norm(seq):
    return np.array(seq)/sum(seq)

color_select=[
    '#34495e',
    '#9b59b6',
    '#c0392b',
    '#f39c12',
    '#f1c40f',
]

for gene in gene_list:
    of=open("host_interaction_%s.csv"%(gene),"w")
    for n in range(25000):
        of.write("%d"%(n+1))
        for ds in ds_keys:
            of.write(",%d"%(gene_host_positions[gene][ds][n]))
        of.write("\n")
    of.close()

site_list=[]
for gene in gene_list:
    for ds in dsets:
        snd_labels,n_labels=snd.label(np.array(smooth(gene_host_locations[gene][ds]))>20)
        objs=snd.find_objects(snd_labels)
        for obj in objs:
            site_list.append([ds,gene,obj[0].start,obj[0].stop])

of=open("site_list.txt","w")
of.write("virus,host_rna,start_v,stop_v,DEN1,ZAFR,ZBRA,ZILM,ZSIN\n")
for site in site_list:
    of.write("%s,%s,%d,%d"%(site[0],site[1],site[2],site[3]))
    for ds in dsets:
        of.write(",%.2f"%(np.array(smooth(gene_host_locations[site[1]][ds])[site[2]:site[3]]).max()))
    of.write("\n")

pairwise_ccf={}
host_cons={}
pairwise_host_ccf={}
cons_val_len={}

global_interactors=[]
global_non_interactors=[]

for gene in gene_list:
    if os.path.exists("phastCons_%s.csv"%(gene)):
        f=open("phastCons_%s.csv"%(gene)).readlines()
        cons_vals=[]
        for i in f:
            l=i.split(",")
            cons_vals.append(float(l[1]))
        host_cons[gene]=[]
        for ds in ds_keys:
            interactions=smooth(gene_host_positions[gene][ds])>5
            interaction_sites=[]
            non_interaction_sites=[]
            for k in range(len(cons_vals)):
                if interactions[k]:
                    interaction_sites.append(cons_vals[k])
                    global_interactors.append(cons_vals[k])
                else:
                    non_interaction_sites.append(cons_vals[k])
                    global_non_interactors.append(cons_vals[k])
            if len(interaction_sites) == 0:
                interaction_sites=[0.0]
            if len(non_interaction_sites) == 0:
                interaction_sites=[0.0]
            host_cons[gene].append(interaction_sites)
            host_cons[gene].append(non_interaction_sites)

        cons_val_len[gene]=len(cons_vals)
        pairwise_host=np.ones((len(ds_keys),len(ds_keys)))
        for i in range(len(ds_keys)):
            for j in range(i+1,len(ds_keys)):
                ccf=np.corrcoef(smooth(gene_host_positions[gene][ds_keys[i]][:len(cons_vals)]),smooth(gene_host_positions[gene][ds_keys[j]][:len(cons_vals)]))[0,1]
                pairwise_host[i,j]=ccf
                pairwise_host[j,i]=ccf
        pairwise_host_ccf[gene]=pairwise_host

    else:
        cons_vals=[]
        print("No info on %s"%(gene))
    pairwise=np.ones((len(ds_keys),len(ds_keys)))
    for i in range(len(ds_keys)):
        for j in range(i+1,len(ds_keys)):
            ccf=np.corrcoef(smooth(gene_host_locations[gene][ds_keys[i]]),smooth(gene_host_locations[gene][ds_keys[j]]))[0,1]
            pairwise[i,j]=ccf
            pairwise[j,i]=ccf
    pairwise_ccf[gene]=pairwise


vp_labels=[]
for ds in ds_keys:
    vp_labels.append("%s\nhost-virus"%(ds))
    vp_labels.append("%s\nno interaction"%(ds))

t_stat=[]
p_stat=[]

for gene in gene_list:
    plt.figure(1,figsize=(30,6))
    plt.subplot2grid((1,13),(0,0),colspan=10)
    plt.title("%s (window smoothing %d nt)"%(gene,window))
    p=[]
    for n,ds in enumerate(dsets):
        p.append(plt.plot(range(11000),smooth(gene_host_locations[gene][ds]),color=color_select[n]))
    plt.xlim(0,11000)
    plt.xticks(range(0,11001,500))
    plt.xlabel("Location (nt)")
    plt.ylabel("SPLASH Read Count")
    plt.legend([k[0] for k in p],dsets.keys())
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.subplot2grid((1,13),(0,11),colspan=2)
    plt.pcolormesh(pairwise_ccf[gene],vmin=-1.0,vmax=1.0,cmap=plt.cm.seismic)
    plt.xticks(np.linspace(0.5,4.5,5),dsets.keys(),rotation='vertical')
    plt.yticks(np.linspace(0.5,4.5,5),dsets.keys())
    plt.colorbar(label='Pearson r')
    plt.savefig("locations_%s.pdf"%(gene),format='pdf',dpi=300)
    plt.clf()

    if gene in host_cons:
        plt.figure(2,figsize=(12,12))
        plt.violinplot(host_cons[gene])
        plt.xticks(range(1,1+len(ds_keys)*2),vp_labels,rotation='vertical')
        pvals=[]
        for i in range(len(ds_keys)):
            t,p=ss.ranksums(host_cons[gene][2*i],host_cons[gene][2*i+1])
            t_stat.append(t)
            p_stat.append(p)
            pvals.append([t,p])
        for i in range(len(pvals)):
            plt.plot([2*i+1+0.1,2*i+1+0.9],[0.81,0.81],color='k')
            plt.annotate("p=%.3g\n(t=%.3f)"%(pvals[i][1],pvals[i][0]),xy=(2*i+1+0.5,0.85),horizontalalignment='center')
        plt.ylabel("PhastCons\n100 Vertebrae")
        plt.title(gene)
        plt.ylim(0,1)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.savefig("conservation_%s.pdf"%(gene),format='pdf',dpi=300)
        plt.clf()

        plt.figure(3,figsize=(30,6))
        plt.subplot2grid((1,13),(0,0),colspan=10)
        plt.title("%s (window smoothing %d nt)"%(gene,window))
        p=[]
        for n,ds in enumerate(dsets):
            p.append(plt.plot(range(cons_val_len[gene]),smooth(gene_host_positions[gene][ds])[:cons_val_len[gene]],color=color_select[n]))
        plt.xlim(0,cons_val_len[gene])
        plt.xticks(range(0,cons_val_len[gene],200))
        plt.xlabel("Location (nt)")
        plt.ylabel("SPLASH Read Count")
        plt.legend([k[0] for k in p],dsets.keys())
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)
        plt.subplot2grid((1,13),(0,11),colspan=2)
        plt.pcolormesh(pairwise_host_ccf[gene],vmin=-1.0,vmax=1.0,cmap=plt.cm.seismic)
        plt.xticks(np.linspace(0.5,4.5,5),dsets.keys(),rotation='vertical')
        plt.yticks(np.linspace(0.5,4.5,5),dsets.keys())
        plt.colorbar(label='Pearson r')
        plt.savefig("host_locations_%s.pdf"%(gene),format='pdf',dpi=300)
        plt.clf()



plt.figure(4,figsize=(30,30))
for n,gene in enumerate(gene_list):
    plt.subplot2grid((14,11),(n//10,n%10))
    plt.title(gene,fontsize='small')
    plt.pcolormesh(pairwise_ccf[gene],vmin=-1.0,vmax=1.0,cmap=plt.cm.seismic)
    if n%10==0:
        plt.yticks(np.linspace(0.5,4.5,5),dsets.keys())
    else:
        plt.yticks(np.linspace(0.5,4.5,5),[])
    if n//10==13:
        plt.xticks(np.linspace(0.5,4.5,5),dsets.keys(),rotation='vertical')
    else:
        plt.xticks(np.linspace(0.5,4.5,5),[])
plt.subplot2grid((14,11),(0,10),rowspan=14)
plt.gca().set_visible(False)
plt.colorbar(ScalarMappable(norm=matplotlib.colors.Normalize(vmin=-1.0,vmax=1.0), cmap=plt.cm.seismic),label="Pearson r")
plt.savefig("all_ccfs.pdf",format='pdf',dpi=300)

plt.figure(5,figsize=(30,30))
for n,gene in enumerate(gene_list):
    if gene in pairwise_host_ccf:
        plt.subplot2grid((14,11),(n//10,n%10))
        plt.title(gene,fontsize='small')
        plt.pcolormesh(pairwise_host_ccf[gene],vmin=-1.0,vmax=1.0,cmap=plt.cm.seismic)
        if n%10==0:
            plt.yticks(np.linspace(0.5,4.5,5),dsets.keys())
        else:
            plt.yticks(np.linspace(0.5,4.5,5),[])
        if n//10==13:
            plt.xticks(np.linspace(0.5,4.5,5),dsets.keys(),rotation='vertical')
        else:
            plt.xticks(np.linspace(0.5,4.5,5),[])
plt.subplot2grid((14,11),(0,10),rowspan=14)
plt.gca().set_visible(False)
plt.colorbar(ScalarMappable(norm=matplotlib.colors.Normalize(vmin=-1.0,vmax=1.0), cmap=plt.cm.seismic),label="Pearson r")
plt.savefig("all_host_ccfs.pdf",format='pdf',dpi=300) 

global_corr_corr={}
global_host_corr={}
denzik_corr_corr={'DEN':[],'ZIK':[]}
denzik_host_corr={'DEN':[],'ZIK':[]}
pairlist=[]
for i in range(len(ds_keys)):
    for j in range(i+1,len(ds_keys)):
        pair_label="%s-%s"%(ds_keys[i],ds_keys[j])
        pairlist.append(pair_label)
        global_corr_corr[pair_label]=[]
        global_host_corr[pair_label]=[]

for n,gene in enumerate(gene_list):
    for i in range(len(ds_keys)):
        for j in range(i+1,len(ds_keys)):
            pair_label="%s-%s"%(ds_keys[i],ds_keys[j])
            global_corr_corr[pair_label].append(pairwise_ccf[gene][i,j])
            if 'DEN' in pair_label:
                denzik_corr_corr['DEN'].append(pairwise_ccf[gene][i,j])
            else:
                denzik_corr_corr['ZIK'].append(pairwise_ccf[gene][i,j])


plt.figure(6,figsize=(25,10))
for n,pl in enumerate(pairlist):
    y,x=np.histogram(global_corr_corr[pl],bins=12,range=(-0.2,1.0))
    x_labels=[]
    for i in range(len(x)-1):
        x_labels.append("%.2f-%.2f"%(x[i],x[i+1]))
    plt.subplot2grid((2,5),(n//5,n%5))
    plt.bar(range(len(y)),y/sum(y),width=0.8,color='k')
    plt.title(pl)
    plt.ylim(0.0,0.6)
    if n//5==1:
        plt.xticks(range(len(y)),x_labels,rotation='vertical')
        plt.xlabel("Pearson r")
    else:
        plt.xticks(range(len(y)),[])
    if n%5==0:
        plt.ylabel("Frequency")
plt.savefig("GlobalVirusCorrelation.pdf",dpi=300,format='pdf')

plt.figure(7,figsize=(24,12))
for n,pl in enumerate(['DEN','ZIK']):
    y,x=np.histogram(denzik_corr_corr[pl],bins=12,range=(-0.2,1.0))
    x_labels=[]
    for i in range(len(x)-1):
        x_labels.append("%.2f-%.2f"%(x[i],x[i+1]))
    plt.subplot2grid((1,2),(0,n))
    plt.bar(range(len(y)),y/sum(y),width=0.8,color='k')
    plt.title(pl)
    plt.ylim(0.0,0.6)
    plt.xticks(range(len(y)),x_labels,rotation='vertical')
    plt.xlabel("Pearson r")
    plt.ylabel("Frequency")
plt.savefig("GlobalDenZikVirusCorrelation.pdf",dpi=300,format='pdf')

for n,gene in enumerate(gene_list):
    if gene in pairwise_host_ccf:
        for i in range(len(ds_keys)):
            for j in range(i+1,len(ds_keys)):
                pair_label="%s-%s"%(ds_keys[i],ds_keys[j])
                global_host_corr[pair_label].append(pairwise_host_ccf[gene][i,j])
                if 'DEN' in pair_label:
                    denzik_host_corr['DEN'].append(pairwise_host_ccf[gene][i,j])
                else:
                    denzik_host_corr['ZIK'].append(pairwise_host_ccf[gene][i,j])

    
plt.figure(8,figsize=(25,10))
for n,pl in enumerate(pairlist):
    y,x=np.histogram(global_host_corr[pl],bins=12,range=(-0.2,1.0))
    x_labels=[]
    for i in range(len(x)-1):
        x_labels.append("%.2f-%.2f"%(x[i],x[i+1]))
    plt.subplot2grid((2,5),(n//5,n%5))
    plt.bar(range(len(y)),y/sum(y),width=0.8,color='k')
    plt.title(pl)
    plt.ylim(0.0,0.6)
    if n//5==1:
        plt.xticks(range(len(y)),x_labels,rotation='vertical')
        plt.xlabel("Pearson r")
    else:
        plt.xticks(range(len(y)),[])
    if n%5==0:
        plt.ylabel("Frequency")
plt.savefig("GlobalHostCorrelation.pdf",dpi=300,format='pdf')
    
plt.figure(9,figsize=(24,12))
for n,pl in enumerate(['DEN','ZIK']):
    y,x=np.histogram(denzik_host_corr[pl],bins=12,range=(-0.2,1.0))
    x_labels=[]
    for i in range(len(x)-1):
        x_labels.append("%.2f-%.2f"%(x[i],x[i+1]))
    plt.subplot2grid((1,2),(0,n))
    plt.bar(range(len(y)),y/sum(y),width=0.8,color='k')
    plt.title(pl)
    plt.ylim(0.0,0.6)
    plt.xticks(range(len(y)),x_labels,rotation='vertical')
    plt.xlabel("Pearson r")
    plt.ylabel("Frequency")
plt.savefig("GlobalDenZikHostCorrelation.pdf",dpi=300,format='pdf')    
    

    





of=open("t_stat.csv","w")
for val in t_stat:
    of.write("%g\n"%(val))
of.close() 

of=open("p_stat.csv","w")
for val in p_stat:
    of.write("%g\n"%(val))
of.close() 


t,p=ss.ranksums(global_interactors,global_non_interactors)

y_inter,x_inter=np.histogram(global_interactors,bins=20,range=(0,1))
y_non_inter,x_non_inter=np.histogram(global_non_interactors,bins=20,range=(0,1))
x_labels=[]
for i in range(len(x_inter)-1):
    x_labels.append("%.2f-%.2f"%(x_inter[i],x_inter[i+1]))

plt.figure(10,figsize=(16,10))

plt.title("t=%.3g, p=%.3g)"%(t,p))
pleg=[]
pleg.append(plt.bar(np.linspace(1,20,20)-0.21,y_inter/sum(y_inter),width=0.4,color='r'))
pleg.append(plt.bar(np.linspace(1,20,20)+0.21,y_non_inter/sum(y_non_inter),width=0.4,color='k'))
plt.legend([p[0] for p in pleg],['Interacting','Non-interacting'])
plt.xticks(np.linspace(1,20,20),x_labels,rotation='vertical')
plt.xlabel("PhastCons 100 Vertebrae")
plt.ylabel("Frequency")

plt.savefig("GlobalConsStatistics.pdf",format='pdf',dpi=300)