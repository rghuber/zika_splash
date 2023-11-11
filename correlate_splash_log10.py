#!/usr/bin/env python3
import sys
import numpy as np
import scipy.stats as ss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
import uncertainties as unc
from scipy.optimize import curve_fit


fpkm={
        "cytoplasm":{},
        "nucleoplasm":{},
        "chromatin":{},
        "full_lysate":{},
        "gene_id":{},
        }

f=open("InfectedGeneAbundance.csv").readlines()
for i in f[1:]:
    l=i.split(",")
    gene=l[1]
    gene_id=l[0]
    t_cytoplasm=float(l[2])
    t_nucleoplasm=float(l[3])
    t_chromatin=float(l[4])
    t_full_lysate=float(l[5])
    if all(np.array([t_cytoplasm,t_nucleoplasm,t_chromatin,t_full_lysate]) > 0):
        fpkm["gene_id"][gene]=gene_id
        fpkm["cytoplasm"][gene]=np.log10(t_cytoplasm)
        fpkm["nucleoplasm"][gene]=np.log10(t_nucleoplasm)
        fpkm["chromatin"][gene]=np.log10(t_chromatin)
        fpkm["full_lysate"][gene]=np.log10(t_full_lysate)

print(fpkm["cytoplasm"]["ERGIC3"])
print(fpkm["nucleoplasm"]["ERGIC3"])
print(fpkm["chromatin"]["ERGIC3"])
print(fpkm["full_lysate"]["ERGIC3"])


abundance={
        "DEN1":{},
        "ZAFR":{},
        "ZBRA":{},
        "ZILM":{},
        "ZSIN":{},
        }

host_peak={
        "DEN1":{},
        "ZAFR":{},
        "ZBRA":{},
        "ZILM":{},
        "ZSIN":{},
        }

denzik_peak={
        "DEN1":{},
        "ZAFR":{},
        "ZBRA":{},
        "ZILM":{},
        "ZSIN":{},
        }

denzik_locations={
        "DEN1":{},
        "ZAFR":{},
        "ZBRA":{},
        "ZILM":{},
        "ZSIN":{},
        }

host_locations={
        "DEN1":{},
        "ZAFR":{},
        "ZBRA":{},
        "ZILM":{},
        "ZSIN":{},
        }

gene_type={
        "DEN1":[],
        "ZAFR":[],
        "ZBRA":[],
        "ZILM":[],
        "ZSIN":[],
}

geneLabel={
        "DEN1":[],
        "ZAFR":[],
        "ZBRA":[],
        "ZILM":[],
        "ZSIN":[],
        }



fpkm_genes=list(fpkm["cytoplasm"].keys())

for ds in abundance:
    f=open("peak_locations_%s.csv"%(ds)).readlines()
    print(ds)
    for i in f[1:]:
        l=i.strip().split(",")
        gene=l[0]
        val=int(l[1])
        dz_peak=int(l[2])
        dz_loc=int(l[3])+1
        h_peak=int(l[4])
        h_loc=int(l[5])+1
        g_type=l[6]

        if (dz_peak >=2) and gene in fpkm_genes:
            abundance[ds][gene]=np.log10(val)
            denzik_peak[ds][gene]=np.log10(dz_peak)
            host_peak[ds][gene]=h_peak
            denzik_locations[ds][gene]=dz_loc
            host_locations[ds][gene]=h_loc
            geneLabel[ds].append(gene)
            gene_type[ds].append(g_type)

print(len(abundance["DEN1"]))
print(len(abundance["ZAFR"]))
print(len(abundance["ZBRA"]))
print(len(abundance["ZILM"]))
print(len(abundance["ZSIN"]))

def f(x,a,b):
    return a*x+b

def predband(x,xd,yd,p,func,conf=0.95):
    x=np.array(x)
    xd=np.array(xd)
    yd=np.array(yd)
    alpha = 1.0 - conf
    N=xd.size
    var_n = len(p)
    q=ss.t.ppf(1.0-alpha/2.0, N-var_n)
    se = np.sqrt(1.0/(N-var_n)*np.sum((yd-func(xd,*p))**2))
    sx = (x-xd.mean())**2
    sxd = np.sum((xd-xd.mean())**2)
    yp=func(x,*p)
    dy=q*se*np.sqrt(1.0+(1.0/N)+(sx/sxd))
    lpb,upb=yp-dy,yp+dy
    return lpb,upb

matplotlib.rcParams.update({
    'savefig.bbox':'tight',
    'font.size':26,
    'axes.linewidth':2,
    'xtick.major.width':2,
    'ytick.major.width':2,
    })

def decide_label(pnt_x,pnt_y,px,upb):
    d=1.e6
    for i in range(len(px)):
        if abs(pnt_x-px[i]) < d:
            d=abs(pnt_x-px[i])
            loc=i
        else:
            loc=i
            break
    if pnt_y > upb[loc]:
        return True
    else:
        return False



for nt,time in enumerate(["cytoplasm","nucleoplasm","chromatin","full_lysate"]):
    for nz,zikv in enumerate(list(abundance.keys())):
        sys.stdout.write("%s%s\n"%(time,zikv))
        sys.stdout.flush()

        of_over=open("overrepresented_splash_%s_%s.csv"%(zikv,time),"w")
        of_over.write("gene_name,log10_splash_count,log10_FPKM,splash_count,FPKM,overrepresentation_sigma,overrepresentation_p\n")

        plt.figure(1,figsize=(12,12))
        x=[]
        y=[]
        for i in geneLabel[zikv]:
            x.append(fpkm[time][i])
            y.append(denzik_peak[zikv][i])
            #y.append(abundance[zikv][i])

        popt,pcov=curve_fit(f,x,y)
        a,b = unc.correlated_values(popt,pcov)
        px=np.linspace(-1.0,5.0,1000)
        py=a*px+b
        nom=unp.nominal_values(py)
        std=unp.std_devs(py)
        lpb,upb=predband(px,x,y,popt,f,conf=0.95)

        labels=[]
        sizes=[5 for i in range(len(x))]
        for i in range(len(x)):
            labels.append(decide_label(x[i],y[i],px,upb))
            if labels[i]:
                sizes[i]=25

        of=open("overrepresented_genes_%s_%s.csv"%(zikv,time),"w")
        of.write("gene,abundance,fpkm,denzik_peak,denzik_peak_location,host_peak,host_peak_location,gene_type\n")

        #plt.subplot2grid((5,4),(nz,nt))
        plt.title("%s-%s r=%.3f"%(zikv,time,np.corrcoef(x,y)[0][1]))
        plt.scatter(x,y,s=sizes,color='k')
        p0=plt.plot(px,nom,color='r',ls='--')
        p1=plt.plot(px,nom-1.96*std,color='grey')
        plt.plot(px,nom+1.96*std,color='grey')
        p2=plt.plot(px,lpb,ls='--',color='orange')
        plt.plot(px,upb,ls='--',color='orange')
        plt.legend([p0[0],p1[0],p2[0]],["Best Fit","95% Confidence Interval","95% Prediction Band"],fontsize='x-small',loc=4)
        for i in range(len(labels)):
            if labels[i]:
                of.write("%s,%d,%f,%d,%d,%d,%d,%s\n"%(
                    geneLabel[zikv][i],
                    10**abundance[zikv][geneLabel[zikv][i]],
                    10**x[i],
                    10**denzik_peak[zikv][geneLabel[zikv][i]],
                    denzik_locations[zikv][geneLabel[zikv][i]],
                    host_peak[zikv][geneLabel[zikv][i]],
                    host_locations[zikv][geneLabel[zikv][i]],
                    gene_type[zikv][i],
                    ))

                dpx=abs(px-x[i])
                loc=np.argmin(dpx)
                dy=y[i]-nom[loc]
                dupb=(upb[loc]-nom[loc])/2
                sigma=dy/dupb

                print("%s %s %s delta=%.3f sigma, p=%.5g"%(zikv,time,geneLabel[zikv][i],dy/dupb,1.0-ss.norm.cdf(sigma,0,1)))
                of_over.write("%s,%f,%f,%d,%f,%f,%g\n"%(geneLabel[zikv][i],y[i],x[i],10**y[i],10**x[i],sigma,1.0-ss.norm.cdf(sigma,0,1)))

                #print(px[loc],x[i],nom[loc],y[i],upb[loc])
                    
                plt.annotate(geneLabel[zikv][i],xy=(x[i],y[i]),fontsize=8)
        plt.xlabel(r"Log$_{10}$ FPKM")
        plt.ylabel(r"Log$_{10}$ SPLASH Peak on %s"%(ds))
        plt.xlim(-1,5)
        plt.ylim(0.0,4.2)
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)

#plt.savefig("ZIKV_SPLASH_vs_FPKM_log10.svg",format='svg',dpi=300)
        plt.savefig("ZIKV_SPLASH_vs_FPKM_log10_%s_%s.pdf"%(zikv,time),format='pdf',dpi=300)
        plt.savefig("ZIKV_SPLASH_vs_FPKM_log10_%s_%s.png"%(zikv,time),format='png',dpi=300)
        plt.clf()
