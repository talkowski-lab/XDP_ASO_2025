#!/apps/lab/miket/anaconda/4.0.5/envs/dg520_py2/bin/python

from __future__ import division
import math
import numpy as np
import shlex, subprocess
import commands
import os, sys


def psl(a,b,c):
    line=[]
    t=a[1:].split(":")
    loc=b+"/"+t[0]+"/"+t[1]+"/"+c+"/trinity_TAF1.psl"
    name=t[2][:-1]
    d0=open(loc,"r").readlines()
    for x in d0:
        tmp=x.split("\t")
        lng=int(tmp[17])
        match=int(tmp[0])
        mismatch=int(tmp[1])
        ngap=int(tmp[4])
        lgap=int(tmp[5])
        if (name == tmp[9]) and (lng>=2) and ((match/(match+mismatch))>0.95) and ngap==0:
            line.append(t[0]+":"+t[1]+":"+name+"\t"+x)
            break

    return line
    
    
def utr(a):           #This function is no longer used in the stage. Transcript will NOT be merged before they can pass junction validation.  
    tr=[]
    process=[]
    m=0
    for x in a:
        n=0
        #if (m in process) ==False:
        trans=x
        for y in a:
            t=trans.split("\t")
            tl1=int(t[1])
            jun1=t[-1].split(",")[1:-1]
            lng1=t[-3].split(",")[1:-2]
            st1=t[9]
            tt=y.split("\t")
            tl2=int(tt[1])
            jun2=tt[-1].split(",")[1:-1]
            lng2=tt[-3].split(",")[1:-2]
            st2=tt[9]
            
            #if t[0]=="capture:Neuron_XDP:TRINITY_DN84_c0_g1_i1" and tt[0]=="rnaseq:NSC_Control:TRINITY_DN37_c1_g1_i2":
            #    print jun1
            #    print jun2
            #    print lng1
            #    print lng2
            #    print dadi
            
            if jun2==jun1 and lng2==lng1 and st1==st2:
                process.append(n)
                if tl2>tl1:
                    trans=y
            #elif (",".join(jun1) in ",".join(jun2)) and ((",".join(lng1) in ",".join(lng2))) and st1==st2:
            #    process.append(n)
            #    trans=y
                    
            n+=1
                
        if (trans in tr) ==False:
            tr.append(trans)
        m+=1
    
    return tr

    
def norm(mat,outf):
    sf=open("/data/talkowski/Samples/XDP/Brain_RNASeq/analysis/DESEq_CAU_OCC/results/sizefactors_matrix_CAU.txt","r").readlines()
    ge=[]
    for x in sf:
        t=x[:-1].split("\t")
        ge.append(float(t[1]))

    
    title=mat[0]
    normfile=open(outf,"w")
    #normfile.write(title)
    for x in mat[1:]:
        t=x[:-1].split("\t")
        val=[0 if y=="NA" else y for y in t[2:]]
        val=map(float,val)
        i=0
        normfile.write(t[0]+"\t"+t[1])
        for y in val:
            normfile.write("\t")
            nv=y/ge[i]
            i+=1
            normfile.write(str(nv))
        normfile.write("\n")
    
    normfile.close()    

    cmd="sort -k1,1 -k2,2 -o "+outf+" "+outf
    stat,out=commands.getstatusoutput(cmd)
    
    noheader=open(outf,"r").readlines()
    normfile=open(outf,"w")
    normfile.writelines([title]+noheader)
    normfile.close() 
    
    return 1


def sj_fpkm(fa,root,sj,ps,phe,cap,cf,cut):
    sjl=open(sj,"r").readlines()
    anno=open(phe,"r").readlines()
    anno1=open(cap,"r").readlines()
    cf1="/".join(cf.split("/")[:-1])+"/pre_alltrans_appearance.txt"
    cf2="/".join(cf.split("/")[:-1])+"/pre_trans_sample.tab"
    cf_counts="/".join(cf.split("/")[:-1])+"/pre_allSamples_SJ.counts.matrix"
    cf_norm="/".join(cf.split("/")[:-1])+"/pre_allSamples_SJ.normalized.matrix"
    rank=phe
    
    
    print "Building transcript junctions..."
    
    w={}
    for x in ps:
        t=x.split("\t")
        w[t[0]]={}
        lg=map(int,t[-3].split(",")[:-1])
        start=map(int,t[-1].split(",")[:-1])
        st=t[9]
        #if set(lg).issuperset(set([121,128])):
        if "121,128" in ",".join(map(str,lg)):
            pos=0
            for y in lg:
                if y==121 and lg[pos+1]==128:
                    anch=pos
                    break
                pos+=1
            lg[anch]=117
            lg[anch+1]=126
            lg=lg[:(anch+1)]+[6]+lg[(anch+1):]
            start[anch+1]=70680717                             #0-based [a,a+length)
            start=start[:(anch+1)]+[70677483]+start[(anch+1):]
        #p=[]
        i=0
        while i<len(start)-1:
            left=start[i]+lg[i]+1   #1-based [a,b]
            right=start[i+1]        #1-based [a,b]
            w[t[0]][str(left)+"-"+str(right)+":"+st]={}

            i+=1
            for x in sjl:
	            #keys=x.split("/")[-2]
                keys0=x.split("/")[-2]
                keys1=x.split("/")[5]
                keys=keys1+"_"+keys0
                w[t[0]][str(left)+"-"+str(right)+":"+st][keys]=0
               
    print "Grepping junction counts from: "
    file_ord=1
    for x in sjl:
        #keys=x.split("/")[-2]
        keys0=x.split("/")[-2]
        keys1=x.split("/")[5]
        keys=keys1+"_"+keys0
        print file_ord,x[:-1].split("/")[-1]
        file_ord+=1
        for y in open(x[:-1],"r").readlines():
            t=y.split("\t")
            if t[3]=="1":
                st="+"
                crd=(t[1])+"-"+(t[2])+":"+st
            elif t[3]=="2":
                st="-"
                crd=(t[1])+"-"+(t[2])+":"+st
            elif t[3]=="0":
                st="."
                crd1=(t[1])+"-"+(t[2])+":+"
                crd2=(t[1])+"-"+(t[2])+":-"
            
            hit=0
            for z in w:
                if st!=".":
                    if t[0]=="X" and (crd in w[z]):
                        w[z][crd][keys]=int(t[6])
                else:
                    if t[0]=="X" and (crd1 in w[z]):
                        w[z][crd1][keys]=int(t[6])
                    elif t[0]=="X" and (crd2 in w[z]):
                        w[z][crd2][keys]=int(t[6])
    
    
    
    print "Building sample-category relationships..."
    
    dic={}
    tot={}
    #dd={}
    ctrl=0
    for x in anno[1:]:
        t=x.split("\t")
        cell=t[1]
        geno=t[2]
	fam="rnaseq_"+cell+"_"+geno
	sa=t[0]
	libsz=int(t[3])
        for y in sjl:
        	if sa == y.split("/")[-2] and ("_RNASeq" in y):
                    sample=y[:-1]
                    if fam in dic:
                        dic[fam].append(sample)
                    else:
                        dic[fam]=[sample]
                    if (sample in tot) == False:
                        tot[sample]=libsz
                    ctrl+=1
                    
    for x in anno1[1:]:
        t=x.split("\t")
        cell=t[1]
        #if cell=="Neurons":
        #    cell="Neuron"
        #elif cell=="Fibroblasts":
        #    cell="Fibro"
        geno=t[2]
        libsz=int(t[3])
        fam="capture_"+cell+"_"+geno
        sa=t[0]
        #libsz=int(t[4][:-1])
        #tot[batch+"_"+sa]=libsz
        for y in sjl:
            if sa == y.split("/")[-2] and ("_CapSeq" in y):
                sample=y[:-1]
                if fam in dic:
                    dic[fam].append(sample)
                else:
                    dic[fam]=[sample]
                if (sample in tot) == False:
                    tot[sample]=libsz
                ctrl+=1
    
    
    
    print "Total", ctrl, "SJ files indexed."
    print "Total", len(tot), "library sizes grepped."         
    #print dadi
    
    
    print "Validating transcripts by junctions..."
    
    winner=[]
    f=open(cf,"w")
    f1=open(cf1,"w")
    f2=open(cf2,"w")
    
    for x in w:          # x: trancript
        app=""
        internal=""
        transa=[]
        transb=[]
        allok=0
        for fam in dic:
            
            fama=[]          # array: each line is a list for junctions counts along a transcript in a certain sample
            famok=0
            for sample in dic[fam]:
                #logf=sample.split("/")[-2]
                logf0=sample.split("/")[-2]
                logf1=sample.split("/")[5]
                logf=logf1+"_"+logf0
                total=tot[sample]
                
                saj=[]
                junc=[]          # list for junction names in their parsing order.
                for y in w[x]:   # y: junction
                    #print w[x][y]
                    #print w[x]
                    saj.append((1000000*w[x][y][logf]/total)+0.00001)
                    junc.append(y)
                    
                if np.min(saj)>float(cut) and ("_RNASeq" in logf):   # only transripts validated in RNASeq counts, not Capture
                #if np.min(saj)>float(cut):
                    famok+=1
                    allok+=1
                    
                    ########
                    parts=fam.split("_")
                    internal+=parts[1]+"_"+parts[2]+"_"+logf+","
                    #########
                    
                if fama==[]:
                    fama=saj
                else:
                    fama=np.vstack([fama,saj])
                    
            #calculate mean, se, num
            if isinstance(fama,list)==False: # normal situation
                trans_avg = np.array(fama).mean(0)
                trans_obs = len(fama)
                trans_sem = np.array(fama).std(0)/np.sqrt(trans_obs)
            else: # if array has only 1 row: only one sample in the category ==> use list instead of array to avoid invalid calculatin of column stats
                trans_avg = fama
                trans_obs = 1
                trans_sem = [0]*len(fama)
            
            #if isinstance(trans_avg,float): #if the transcript only have 1 junction: 2 exons ==> force column stats result into a list format
            #    trans_avg=[trans_avg]
            #    trans_sem=[trans_sem]

            ii=0
            while ii < len(junc):
                content=x+"\t"+fam+"\t"+junc[ii]+"\t"+str(trans_avg[ii])+"\t"+str(trans_sem[ii])+"\t"+str(trans_obs)+"\n"
                transa.append(content) 
                ii+=1
            
            if "70644089-70659646:+" in w[x]:
                print fam,famok,len(dic[fam])
            
            #if famok >= 1:
            if famok >= 0.5*len(dic[fam]):
                app+=fam+","
            
        #write to junction counts output of a certain transcript for R
        if internal!="":
            f2.write(x+"\t"+internal[:-1]+"\n")    
        if app!="":
            f.writelines(transa)
            f1.write(x+"\t"+app[:-1]+"\n")      
            winner.append(x)      
    
    f.close()
    f1.close()
    f2.close()
    
    dd=open(rank,"r").readlines()
    dd=dd[1:]
    
    countlist=[]
    print winner
    f1=open(cf_counts,"w")
    for x in winner:
        for y in w[x]:
            header="Transcript\tJunction"
            line=x+"\t"+y
            for z in dd:
                #if "rnaseq" in z:
                t=z.split("\t")[0]
                tmp="Brain_RNASeq_"+t
                line+="\t"+str(w[x][y][tmp])
                header+="\t"+t
            countlist.append(line+"\n")
            
    
    countlist=[header+"\n"]+sorted(countlist)
    f1.writelines(countlist)
    f1.close()
    
    #cmd="sort -k1,1 -k2,2 -o "+cf_counts+" "+cf_counts
    #stat,out=commands.getstatusoutput(cmd)
    
    norm(countlist,cf_norm)

    return winner
    
 
def cmp(a,b,inp1,inp2,inp3,fp,cutoff):
    p=[]
    q=[]
    folder="/".join(fp.split("/")[:-1])
    approach=fp.split("/")[-3]
    approach1=a.split("/")[-2]
    #if approach1!=approach:
    #    print approach1, "does not match", approach
    outpFA=folder+"/pre_TAF1_contigs.fa"
    outpPS=folder+"/pre_TAF1_contigs.psl"
    d=open(a,"r").readlines()
    i=0
    while (i+1)<len(d):
        j=psl(d[i],b,approach)
        if j!=[]:
            p=p+j
            
        i+=2
    
    p1=p
    #p1=utr(p)
    #print len(p1),"transcripts left after merging."
    #for xx in p1:
    #    print xx
    #print dadi

    p3=sj_fpkm(a,b,inp1,p1,inp2,inp3,fp,cutoff)
    print len(set(p3)),"transcripts passed validation."
    print "Preparing the final FASTA and PSL files..."
    
    i=0
    while (i+1)<len(d):
        fai=d[i][1:-1]
        for x in p3:
            ids=x.split("\t")[0]
            if ids==fai:
                q.append(d[i])
                q.append(d[i+1])
                break
            
        i+=2   
    
    f=open(outpPS,"w")    
    for x in p1:
        t=x.split("\t")
        lg=map(int,t[-3].split(",")[:-1])
        start=map(int,t[-1].split(",")[:-1])
        if t[0] in p3:
            #if set(lg).issuperset(set([121,128])):
            if "121,128" in ",".join(map(str,lg)):
                pos=0
                for y in lg:
                    if y==121 and lg[pos+1]==128:
                        anch=pos
                        break
                    pos+=1
                lg[anch]=117
                lg[anch+1]=126
                lg=lg[:(anch+1)]+[6]+lg[(anch+1):]
                start[anch+1]=70680717
                start=start[:(anch+1)]+[70677483]+start[(anch+1):] 
                lg=map(str,lg)
                start=map(str,start)
                t[-3]=",".join(lg)+","
                t[-1]=",".join(start)+","
                f.write("\t".join(t)+"\n")
            else:
                f.writelines(x)
    f.close()
    
    f=open(outpFA,"w")
    f.writelines(q)
    f.close()
        
    return 1
        
cmp("/data/talkowski/Samples/XDP/Brain_RNASeq/data/Trinity/RNA_CAP_mergebam/trinity_CAU_mergedbams_trimmed.dedup.bam/Trinity_1line.fasta",
    "/data/talkowski/Samples/XDP/Brain_RNASeq/data/Trinity/RNA_CAP_mergebam",
    "/data/talkowski/Samples/XDP/Brain_RNASeq/data/Trinity/RNA_CAP_merge/CAU_SJ_list",
    "/data/talkowski/Samples/XDP/Brain_RNASeq/data/Trinity/RNA_CAP_merge/Design_matrix_CAU_RNASeq.txt",
    "/data/talkowski/Samples/XDP/Brain_RNASeq/data/Trinity/RNA_CAP_merge/Design_matrix_CAU_CapSeq.txt",
    "/data/talkowski/Samples/XDP/Brain_RNASeq/data/Trinity/RNA_CAP_mergebam/Caudate/all/pre_allSamples_SJ.fpkm.processed.matrix",
    "0.00001")

#Final version
#cmp("/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_DN/contigs-3.fa",
#    "/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/",
#    "/data/talkowski/dg520/projects/xdp/output/SJfiles2ndpass_rnaseq_up2Jun17_Uniform_all.txt",
#    "/data/talkowski/dg520/projects/xdp/output/phenoSamples_lessFriendly_finalUsage_all.txt",
#    "/data/talkowski/dg520/projects/xdp/output/captureSample_withCarrier.txt",
#    "/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_DN/cpm00/pre_allSamples_SJ.fpkm.processed.matrix",
#    "0.00001")

#cmp("/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_GG/contigs-3.fa",
#    "/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/",
#    "/data/talkowski/dg520/projects/xdp/output/SJfiles2ndpass_rnaseq_up2Jun17_Uniform_all.txt",
#    "/data/talkowski/dg520/projects/xdp/output/phenoSamples_lessFriendly_finalUsage_all.txt",
#    "/data/talkowski/dg520/projects/xdp/output/captureSample_withCarrier.txt",
#    "/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_GG/cpm00/pre_allSamples_SJ.fpkm.processed.matrix",
#    "0.00001")
    
#cmp("/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_DN_dedup/contigs-3.fa",
#    "/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/",
#    "/data/talkowski/dg520/projects/xdp/output/SJfiles2ndpass_rnaseq_up2Jun17_Uniform_all.txt",
#    "/data/talkowski/dg520/projects/xdp/output/phenoSamples_lessFriendly_finalUsage_all.txt",
#    "/data/talkowski/dg520/projects/xdp/output/captureSample_withCarrier.txt",
#    "/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_DN_dedup/cpm00/pre_allSamples_SJ.fpkm.processed.matrix",
#    "0.00001")

#cmp("/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_GG_dedup/contigs-3.fa",
#    "/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/",
#    "/data/talkowski/dg520/projects/xdp/output/SJfiles2ndpass_rnaseq_up2Jun17_Uniform_all.txt",
#    "/data/talkowski/dg520/projects/xdp/output/phenoSamples_lessFriendly_finalUsage_all.txt",
#    "/data/talkowski/dg520/projects/xdp/output/captureSample_withCarrier.txt",
#    "/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_GG_dedup/cpm00/pre_allSamples_SJ.fpkm.processed.matrix",
#    "0.00001")
    

#bsub -sla miket_sc -q medium -J test "/data/talkowski/dg520/projects/xdp/scripts/validate_transcripts_and_count_junctions.py"
