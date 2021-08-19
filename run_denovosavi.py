#!/usr/bin/env python

import pandas as pd
from scipy.stats import binom
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib as mpl
import math
import seaborn as sns
pd.options.mode.chained_assignment = None



def main():

    # Check usage
    if len(sys.argv) not in [4]:
        sys.exit("Usage: python denovosavi.py file1 file2 file3 file4")

    # Parse command-line arguments
    input1 = sys.argv[1]
    input2 = sys.argv[2]
    input3 = sys.argv[3]
    input4 = sys.argv[4]

    
    #SAVI data read-in
    #read PDfilter information
    inPD2, inPD3 = readsavi(input1, input2)
    #read parent variant information
    inFM2, inFM3 = readsavi(input3, input4)

    #add parent variant info to PD variants
    FM2specific = [x for x in inFM2.columns if x not in inPD2.columns]
    FM3specific = [x for x in inFM3.columns if x not in inPD3.columns]

    mer2 = inPD2.apply(lambda x: addparentinfo(x,inFM2,FM2specific),axis=1)
    mer2 = pd.concat([inPD2, mer2], axis=1)

    mer3 = inPD3.apply(lambda x: addparentinfo(x,inFM3,FM3specific),axis=1)
    mer3 = pd.concat([inPD3, mer3], axis=1)

    #merge mer2 and mer3
    mer23 = pd.concat([mer2,mer3])
    # print(len(mer23.columns),len(mer2.columns),len(mer3.columns))
    mer = mer23[mer3.columns]

    
    
    #Get denovo
    pre_dnm = getdnm(mer)
    print("preliminary n of DNM:", len(pre_dnm))

    #calculate binom p for heterozygous DNMs
    calc_binom(pre_dnm)

    #adjust binom p
    adj = p_adjust_bh(pre_dnm["binom_pval"].tolist())
    pre_dnm.loc[:,"bh-adj_pval"]=adj

    #filter by adjust_p
    dnm = pre_dnm.loc[pre_dnm["bh-adj_pval"] > 0.01,:]
    print("final DNMs:", len(dnm))

    #dnm graph
    dnm_graph(dnm)
    
    
    
    
    #get NH-SNP
    mer["F_freq"] = mer["F_freq"].astype("int")
    nhsnp = hmz(mer)

    #common NH-SNP
    com = common_hmz(nhsnp)

    #non-common NH-SNP
    noncom = noncommon_hmz(nhsnp)




    #calculate ratio of functional class
    mis, syn, ns, sp, other = calc_ratio(com)
    print("missense/synonymous ratio for common NH-SNP = (", mis, "/", syn, ")")
    plot_pie([mis, syn, ns, sp, other],"common_NH_SNP_pie")

    mis, syn, ns, sp, other = calc_ratio(noncom)
    print("missense/synonymous ratio for non-common NH-SNP = (", mis, "/", syn, ")")
    plot_pie([mis, syn, ns, sp, other],"noncommon_NH_SNP_pie")

    mis, syn, ns, sp, other = calc_ratio(dnm)
    print("missense/synonymous ratio for DNM = (", mis, "/", syn, ")")
    plot_pie([mis, syn, ns, sp, other],"DNM_pie")

    
    
    
    
def replace2(col):
    col = [x.replace("_1", "_FM") for x in col]
    col = [x.replace("_2", "_BLO") for x in col]
    col = [x.replace("_3", "_F") for x in col]
    col = [x.replace("_4", "_M") for x in col]
    col = [x.replace("1_", "FM_") for x in col]
    col = [x.replace("2_", "BLO_") for x in col]
    col = [x.replace("3_", "F_") for x in col]
    col = [x.replace("4_", "M_") for x in col]
    col = [x.replace("2-", "BLO-") for x in col]
    col = [x.replace("3-", "F-") for x in col]
    col = [x.replace("4-", "M-") for x in col]
    return col
    
def replace3(col):
    col = [x.replace("_1", "_FM") for x in col]
    col = [x.replace("_2", "_BLO") for x in col]
    col = [x.replace("_3", "_AVM") for x in col]
    col = [x.replace("_4", "_F") for x in col]
    col = [x.replace("_5", "_M") for x in col]
    col = [x.replace("1_", "FM_") for x in col]
    col = [x.replace("2_", "BLO_") for x in col]
    col = [x.replace("3_", "AVM_") for x in col]
    col = [x.replace("4_", "F_") for x in col]
    col = [x.replace("5_", "M_") for x in col]
    col = [x.replace("2-", "BLO-") for x in col]
    col = [x.replace("3-", "AVM-") for x in col]
    col = [x.replace("4-", "F-") for x in col]
    col = [x.replace("5-", "M-") for x in col]
    return col

def readsavi(file1,file2):
    """read output from SAVI report and change column names"""
    
    in2 = pd.read_csv(file1,sep="\t",header=0,low_memory=False)
    in3 = pd.read_csv(file2,sep="\t",header=0,low_memory=False)

    in2 = in2[[x for x in in2.columns if "Unnamed:" not in x]]
    in3 = in3[[x for x in in3.columns if "Unnamed:" not in x]]

    newcol = in2.columns.values
    newcol[0] = "ID"
    in2.columns = newcol

    newcol = in3.columns.values
    newcol[0] = "ID"
    in3.columns = newcol

    in2.columns = replace2(in2.columns.tolist())
    in3.columns = replace3(in3.columns.tolist())

    return in2, in3


def addparentinfo(inPD,inFM,FMspecific):
    """insert parent information in PDfilter variants"""
    
    K = inFM.loc[(inFM["ID"] == (str(inPD["ID"])+"_FM"))&
                  (inFM["#chromosome"] == inPD["#chromosome"])&
              (inFM["position"] == inPD["position"])&
              (inFM["alt"] == inPD["alt"]),FMspecific]
    return K.iloc[0,:]

def getdnm(df):
    """returns potential DNMs"""
    podnm = df.loc[(df["totdepth_F"] >= 10)
                &(df["totdepth_M"] >= 10)
                 &(df["totdepth_BLO"] >= 20)
                &(df["FM_freq"]==0)
                &(df["BLO_freq"]!=0)
                &(df["meganormal_id"]=="-")
                &(df["strand_bias_BLO"] > 0.01)
                , :]
    return podnm



def calc_binom(pre_dnm):
    """whether VAF follows binomial distribution with n(totdepth) and p=0.5 (probability) 
        binom.cdf(alt read, totdepth, 0.5)
        output= pval of whether the alt read is within the totdepth and p=0.5 binomial distribution
        if pval < 0.05 or >0.95 then it is not within the distribution"""

    for i in pre_dnm.index:
        depth = pre_dnm.loc[i,"totdepth_BLO"]
        alt = pre_dnm.loc[i,"altdepth_BLO"]
        if alt <= (depth/2):
            pval=binom.cdf(alt,depth,0.5)
        else:
            pval=1-(binom.cdf(alt-1,depth,0.5))
    
        pre_dnm.loc[i,"binom_pval"] = pval
        

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def hmz(df):
    """return variants that are not homozygous in patient but not in parents(NH-SNP) """ 
    hm = df.loc[
               (df["totdepth_BLO"] >= 20)
              &(df["totdepth_F"] >= 20)
              &(df["totdepth_M"] >= 20)
                
              &(df["BLO_freq"] > 85)
              &(df["F_freq"] <= 85)
              &(df["M_freq"] <= 85)
    
              &(df["meganormal_id"]=="-")
    
    
              &(df["strand_bias_BLO"] > 0.01)
               , :]
    
    #remove chr X
    hm = hm.loc[hm["#chromosome"]!="chrX",:]
    return hm
   
def common_hmz(all_hm):
    """return common NH-SNP """
    com_hz = all_hm.loc[(all_hm["snp.common"].str.contains('1')),:]
    return com_hz

def noncommon_hmz(all_hm):
    """return common NH-SNP """
    noncom_hz = all_hm.loc[(~all_hm["snp.common"].str.contains('1')),:]
    return noncom_hz

def calc_ratio(df):
    """calculate the ratio of functional class of the identified variants"""
    fc = df["Functional_Class"].value_counts()

    mis = 0
    syn = 0
    ns = 0
    sp = 0
    other = 0
    
    for typ in fc.index:
        if ("NONSENSE" in typ):
            ns = ns + fc[typ]
        elif ("NONSENSE" not in typ) and ("MISSENSE" in typ):
            mis = mis+ fc[typ]
        elif ("NONSENSE" not in typ) and ("MISSENSE" not in typ) and ("SILENT" in typ) :
            syn = syn+ fc[typ]
        else:
            typlist = df.loc[df["Functional_Class"] == typ,"Effect"].tolist()
            for i in typlist:
                if "synonymous" in i:
                    syn = syn+1
                elif ("synonymous" not in i) and ("splice" in i):
                    sp = sp+1
                elif ("synonymous" not in i) and ("splice" not in i):
                    other = other + 1
                    #print("variant type included in other: ", i)
        
    if mis+syn+ns+sp+other > len(df):
        print("duplicate detected")
        
    return mis, syn, ns, sp, other

def dnm_graph(dnm):
    """draw histogram showing the number of denovo mutations per patient"""
    if len(dnm) != 0:
        
        pats = dnm["ID"].value_counts()
        
        label_size = 20
        mpl.rcParams['xtick.labelsize'] = label_size 
        mpl.rcParams['ytick.labelsize'] = label_size
        mpl.rcParams['figure.figsize'] = [8,8]

        ax =  plt.subplot(111)

        x = []
        y = []
        
        plt.hist(pats,color="#34675C")
        plt.ylabel("number of patients", size = label_size)
        plt.xlabel("number of denovo mutation/patient", size = label_size)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.savefig("denovoSAVI_fam_hist.pdf")
        plt.show()
    else:
        print("dnm_graph: no variants to plot")

def plot_pie(numbers,data_type):
    """draw pie chart showing types of variants"""
    if sum(numbers) >0 :
    
        mpl.rcParams['figure.figsize'] = [4,4]
        fig, ax =  plt.subplots()
        labels = ["Missense", "Synonymous","Nonsense","Splice region","Other"]
        color = ["#f26d5b","#3b8686","#FFBC42","#537ac9","#c6a49a"]
        patches, texts  = plt.pie(numbers,# labels=labels, 
                          labeldistance=1.1,
                              colors=color,textprops={'fontsize': 5})
    
        #set transparency in color
        for i in patches:
            i.set_alpha(1)

        plt.legend(patches, labels,loc="center left",
              bbox_to_anchor=(1, 0, 0, 1), fontsize=10,frameon=False)
        
        plt.savefig("denovoSAVI_"+data_type+".pdf", bbox_inches="tight")
        plt.show()
    else:
        print("plot_pie: no variants to plot")
        
        
#-------------------------------   
if __name__ == "__main__":
    main()

