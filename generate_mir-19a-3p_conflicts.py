#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mp
import matplotlib.collections as mc
from matplotlib.cm import ScalarMappable
import scipy.ndimage as snd

matplotlib.rcParams.update({
    'savefig.bbox':'tight',
    'font.size':26,
    'axes.linewidth':2,
    'xtick.major.width':2,
    'ytick.major.width':2,
})


def arc_plot(clean_pairs,ax1=None,color='b'):

    patches=[]

    for i in clean_pairs:
        radius=abs(float(i[1])-float(i[0]))/2.0
        center=abs(float(i[1])+float(i[0]))/2.0
        patches.append(mp.Wedge([center,0.0],radius,0.0,180.0,width=10.0))

    collection=mc.PatchCollection(patches,edgecolor='none',facecolor=color,alpha=0.5)
    ax1.add_collection(collection)

def arc_plot_inv(clean_pairs,ax1=None,color='b'):

    patches=[]

    for i in clean_pairs:
        if i[0] > 850 and i[0] < 1050 and i[1] > 850 and i[1] < 1050:
            radius=abs(float(i[1])-float(i[0]))/2.0
            center=abs(float(i[1])+float(i[0]))/2.0
            patches.append(mp.Wedge([center,0.0],radius,180.0,360.0,width=5.0))

    collection=mc.PatchCollection(patches,edgecolor='none',facecolor=color,alpha=0.5)
    ax1.add_collection(collection)

query_location=range(9684//10,9700//10)

dsets={
    "DEN1":"EU081230",
    "ZAFR":"AY632535",
    "ZBRA":"KU497555",
    "ZILM":"KJ776791",
    "ZSIN":"KY241786",
}

ds_keys=list(dsets.keys())

projmaps={n:np.zeros((1100,1100)) for n in dsets}

for ds in dsets:
    vname=dsets[ds]
    f=open("../%s.chimera-flt.csv"%(ds)).readlines()
    for i in f:
        l=i.split()
        gene_a=l[2]
        start_a=(int(l[3])-1)//10
        stop_a=(int(l[4])-1)//10
        gene_b=l[9]
        start_b=(int(l[10])-1)//10
        stop_b=(int(l[11])-1)//10
        if gene_a == vname and gene_b == vname:
            for x in range(start_a,stop_a):
                for y in range(start_b,stop_b):
                    projmaps[ds][x,y]+=1
                    projmaps[ds][y,x]+=1

all_conflict=[]
for ds in dsets:
    print(ds)
    conflict=np.zeros(1100)
    for loc in query_location:
        conflict+=projmaps[ds][loc,:]
    conflict_l,n_labels=snd.label(conflict)
    all_conflict.append(conflict/max(conflict))
    objs=snd.find_objects(conflict_l)
    for obj in objs:
        if max(conflict[obj[0]]) >= 10:
            print("%d-%d:%d"%(obj[0].start*10,obj[0].stop*10,max(conflict[obj[0]])))
    print()

for ds in dsets:
    #sys.stdout.write("%s\n"%(ds))
    for i in range(len(projmaps[ds])):
        for j in range(i,len(projmaps[ds])):
            projmaps[ds][i,j]=0
            #sys.stdout.write(" %3d"%(projmaps[ds][i,j]))
        #sys.stdout.write("\n")


zilm_l,n_labels=snd.label(projmaps['ZILM']>=3)
zilm_objs=snd.find_objects(zilm_l)
clean_pairs=[]
of=open("table_miR-19a-3p.csv","w")
of.write("lhs_start,lhs_stop,rhs_start,rhs_stop,n_SPLASH_reads\n")
for obj in zilm_objs:
    if 968 in range(obj[0].start,obj[0].stop) or 968 in range(obj[1].start,obj[1].stop):
        of.write("%d,%d,%d,%d,%d\n"%(obj[0].start*10,obj[0].stop*10,obj[1].start*10,obj[1].stop*10,projmaps['ZILM'][obj].max()))
    if projmaps['ZILM'][obj].max()>=10:
        clean_pairs.append([np.mean([obj[0].start,obj[0].stop]),np.mean([obj[1].start,obj[1].stop])])
print(clean_pairs)


plt.figure(1,figsize=(30,15))

plt.subplot2grid((4,1),(0,0),rowspan=2)
ax1=plt.gca()
arc_plot(clean_pairs,ax1,'grey')
plt.bar([968],[550],color='k',width=2,zorder=5)
plt.annotate("miR-19a-3p",xy=(978,500))
plt.xlim(0,1100)
plt.ylim(0,550)
plt.xticks(range(0,1101,50),[])
plt.yticks([])
plt.grid(True)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.title("Virus-Virus interactions competing with Virus - mir-19a-3p")

plt.subplot2grid((4,1),(2,0))
plt.pcolormesh(np.array(all_conflict),vmin=0,vmax=1.0,cmap=plt.cm.Reds,zorder=4)
plt.bar([968],[5],color='k',width=2,zorder=5)
plt.xticks(range(0,1101,50),[])
plt.yticks(np.linspace(0.5,4.5,5),dsets.keys())
plt.grid(True)

plt.subplot2grid((4,1),(3,0))
ax1=plt.gca()
arc_plot_inv(clean_pairs,ax1,'r')
plt.bar([968],[-35],color='k',width=2,zorder=5)
plt.xlim(0,1100)
plt.ylim(-35,0)
plt.xticks(range(0,1101,50),range(0,11001,500),rotation='vertical')
plt.yticks([])
plt.grid(True)
plt.gca().spines['bottom'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

plt.savefig("mir_19a-3p.pdf",format='pdf',dpi=300)