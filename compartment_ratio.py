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
    'cytoplasm 1',
    'cytoplasm 2',
    'nucleoplasm 1',
    'nucleoplasm 2',
    'chromatin 1',
    'chromatin 2',
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
    #for n in ['infected','non_infected']:
    #    for m in ['cytoplasm','nucleoplasm','chromatin','full_lysate']:
    #        data_table[gene_id][n][m]=np.mean(data_table[gene_id][n][m])
    
#print(data_table[gene_ids[0]])

print("gene_id,gene_name,inf_mean(log(c/c+n)),inf_std(log(c/c+n)),noninf_mean(log(c/c+n)),noninf_std(log(c/c+n)),t,p")

heatmap=[]
heatmap_labels=[]

cytoplasm_nucleoplasm_ratio={'infected':{},'non-infected':{}}
for gene in gene_ids:
    
    inf_cp=data_table[gene]['infected']['cytoplasm']
    inf_np=data_table[gene]['infected']['nucleoplasm']
    non_cp=data_table[gene]['non_infected']['cytoplasm']
    non_np=data_table[gene]['non_infected']['nucleoplasm']
    
    cytoplasm_nucleoplasm_ratio_infected=np.mean(inf_cp)/(np.mean(inf_cp)+np.mean(inf_np))
    mean_denom_inf=np.mean(inf_cp)+np.mean(inf_np)
    std_denom_inf=(np.std(inf_cp)**2.0+np.std(inf_np)**2)**0.5
    std_inf=abs(cytoplasm_nucleoplasm_ratio_infected)*((np.std(inf_cp)/np.mean(inf_cp))**2.0+(std_denom_inf/mean_denom_inf)**2.0)**0.5


    cytoplasm_nucleoplasm_ratio_noninfected=np.mean(non_cp)/(np.mean(non_cp)+np.mean(non_np))
    mean_denom_non=np.mean(non_cp)+np.mean(non_np)
    std_denom_non=(np.std(non_cp)**2.0+np.std(non_np)**2)**0.5
    std_non=abs(cytoplasm_nucleoplasm_ratio_noninfected)*((np.std(non_cp)/np.mean(non_cp))**2.0+(std_denom_non/mean_denom_non)**2.0)**0.5

    t,p=ss.ttest_ind_from_stats(cytoplasm_nucleoplasm_ratio_infected,std_inf,2,cytoplasm_nucleoplasm_ratio_noninfected,std_non,2,equal_var=False)

    if p < 0.05 and abs(cytoplasm_nucleoplasm_ratio_infected-cytoplasm_nucleoplasm_ratio_noninfected) > 0.2:
        print("%s,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f"%(gene,data_table[gene]['name'],cytoplasm_nucleoplasm_ratio_infected,std_inf,cytoplasm_nucleoplasm_ratio_noninfected,std_non,t,p))

        heatmap.append([cytoplasm_nucleoplasm_ratio_noninfected,cytoplasm_nucleoplasm_ratio_infected])
        if data_table[gene]['name'] != '-':
            heatmap_labels.append(data_table[gene]['name'])
        else:
            heatmap_labels.append(gene)

        #plt.figure(1,figsize=(6,12))
        #plt.title("%s (%s) p=%.3f"%(gene,data_table[gene]['name'],p))
        #b=plt.bar([1,2],[cytoplasm_nucleoplasm_ratio_noninfected,cytoplasm_nucleoplasm_ratio_infected],yerr=[std_non,std_inf],align='center',ecolor='k')
        #b[0].set_color('#28B463')
        #b[1].set_color('#E74C3C')
        #plt.xticks([1,2],['Non-Infected','Infected'])
        #plt.ylabel(r"$\frac{TPM_{Cytoplasm}}{TPM_{Cytoplasm}+TPM_{Nucleoplasm}}$",fontsize='x-large')
        #plt.ylim(0,1)
        #plt.savefig("CompartmentShift/Shift_%s.pdf"%(gene),format='pdf',dpi=300)
        #plt.clf()



heatmap=np.array(heatmap)

order=np.argsort(heatmap[:,1]-heatmap[:,0])

heatmap_sorted=[heatmap[n] for n in order]
heatmap_labels_sorted=[heatmap_labels[n] for n in order]

plt.figure(1,figsize=(8,40))
plt.pcolormesh(heatmap_sorted,cmap=plt.cm.seismic,vmin=0,vmax=1)
plt.yticks(np.linspace(0.5,len(heatmap_labels)-0.5,len(heatmap_labels)),heatmap_labels_sorted,fontsize='x-small')
plt.xticks([0.5,1.5],['Non-Infected','Infected'])
plt.colorbar(label=r"$\frac{TPM_{Cytoplasm}}{TPM_{Cytoplasm}+TPM_{Nucleoplasm}}$")
plt.savefig("CompartmentShift_Heatmap.pdf",format='pdf',dpi=300)