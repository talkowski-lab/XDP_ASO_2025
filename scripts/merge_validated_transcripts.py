#!/apps/lab/miket/anaconda/4.0.5/envs/dg520_py2/bin/python

from __future__ import division
import math
import numpy as np
import shlex, subprocess
import commands
import os, sys



def new(left,right,genepred,psl,root):
    guide_out=root+"/new_and_known_transcripts.txt"
    ref=open(genepred,"r").readlines()
    
    p=[]
    q=[]
    
    for y in psl:
        hit=0
        startl1=psl[y]["jun"]
        lg1=psl[y]["len"]
        endl1=[m+n for m,n in zip(startl1,lg1)]
        startl1=map(str,startl1)
        endl1=map(str,endl1)
        st1=psl[y]["str"]
        for x in ref:
            t=x.split("\t")
            ch=t[1]
            if ch=="X":
                startl=t[8].split(",")[:-1]
                endl=t[9].split(",")[:-1]
                st=t[2]
                if int(startl[0])>left and int(endl[-1])<right:
                    if len(startl1)>1:     
                        v1=",".join(startl[1:])
                        v2=",".join(endl[:-1])
                        v3=",".join(startl1[1:])
                        v4=",".join(endl1[:-1])
                    
                        if v1==v3 and v2==v4 and st==st1:
                            p.append(y+"\t"+t[0])
                            hit=1
                            break
                    elif len(startl1)==1:
                        if len(startl)==1:
                            if max(int(startl1[0]),int(startl[0]))<min(int(endl1[0]),int(endl[0])) and st==st1:
                                p.append(y+"\t"+t[0])
                                hit+=1
                                break
                        
                        
        if hit==0:
            q.append(y+"\t.")
    
    newid={}        
    f=open(guide_out,"w")
    i=0
    for x in p:
        oldname=x.split("\t")[0]
        newname=x.split("\t")[1]
        if newname=="ENST00000276072":
            f.write(x+"\tcTAF1\n")
            newid[oldname]="cTAF1"
        elif newname=="ENST00000423759":
            f.write(x+"\tnTAF1\n")
            newid[oldname]="nTAF1"
        else:
            i+=1
            f.write(x+"\tKnown"+str(i)+"\n")
            newid[oldname]="Known"+str(i)
    i=0
    for x in q:
        i+=1
        f.write(x+"\tNew"+str(i)+"\n")
        oldname=x.split("\t")[0]
        newid[oldname]="New"+str(i)
    f.close()
    
    return newid


def getsj(oldnameList,newnameList,root,inp,outp):
    d=open(root+"/"+inp,"r").readlines()
    f1=open(root+"/"+outp,"w")
    f1.write(d[0])
    for x in d[1:]:
        t=x.split("\t")
        if t[0] in oldnameList:
            f1.write(newnameList[t[0]]+"\t"+"\t".join(t[1:]))
    f1.close()
        
    return 1
    

def judge(a,b):
    judgement=0
    if len(a["jun"])==len(b["jun"]) and a["jun"][1:]==b["jun"][1:] and a["len"][1:-1]==b["len"][1:-1]:  # only UTR difference
        judgement=1
    elif len(a["jun"])!=len(b["jun"]):           #contain or not contain
        if len(a["jun"])>len(b["jun"]):
            long =a
            short=b
        else:
            long =b
            short=a
        j1=map(str,long["jun"][1:])
        j2=map(str,short["jun"][1:])
        j3=map(str,long["len"][1:-1])
        j4=map(str,short["len"][1:-1])
        if ((",".join(j2) in ",".join(j1)) and min(len(a["len"]),len(b["len"]))==2) or ((",".join(j2) in ",".join(j1)) and (",".join(j4) in ",".join(j3)) and min(len(a["len"]),len(b["len"]))>=2):         #long one contains the short one in terms of junctions
            judgement=1
            
            #if short["jun"][1]==long["jun"][1]:
            #    p=long["jun"].index(short["jun"][-1])
            #    q=short["jun"].index(short["jun"][-1])
            #    if short["len"][q] <= long["len"][p]:
            #        judgement=1
            #elif short["jun"][-1]==long["jun"][-1]:
            #    p=long["jun"].index(short["jun"][1])
            #    q=short["jun"].index(short["jun"][1])
            #    if short["len"][q-1] <= long["len"][p-1]:
            #        judgement=1
            #else:
            #    p1=long["jun"].index(short["jun"][1])
            #    q1=short["jun"].index(short["jun"][1])
            #    p2=long["jun"].index(short["jun"][-1])
            #    q2=short["jun"].index(short["jun"][-1])
            #    if (short["len"][q1-1] <= long["len"][p1-1]) and (short["len"][q2] <= long["len"][p2]):
            #        judgement=1
        
    return judgement


def cmp(a,b):
    a["inf"][18] = str(len(a["jun"]))   #correct number of exon first
    a["inf"][7]  = str(len(a["jun"])-1) #correct number of intron first
    b["inf"][18] = str(len(b["jun"]))   #correct number of exon first
    b["inf"][7]  = str(len(b["jun"])-1) #correct number of intron first
        
    precondition=judge(a,b)
    res=0
    if a["chr"]==b["chr"] and a["str"]==b["str"] and precondition==1 and sorted(a["tis"])==sorted(b["tis"]):
    #if a["chr"]==b["chr"] and a["str"]==b["str"] and precondition==1:
        if a["inf"][0]=="TRINITY_GG_13_c215_g2_i2" or b["inf"][0]=="TRINITY_GG_13_c215_g2_i2":
            print a["seq"][0:10],len(a["jun"]),a["jun"][1],a["jun"][-1],a["inf"][0]
            print b["seq"][0:10],len(b["jun"]),b["jun"][1],b["jun"][-1],b["inf"][0]
        res=1
        if len(a["jun"]) < len(b["jun"]):
            base=b
            alt =a
        elif len(a["jun"]) > len(b["jun"]):
            base=a
            alt =b
        else:
            if int(a["totalLen"])>=int(b["totalLen"]):
                base=b
                alt =a
            else:
                base=a
                alt =b
            
        
        #newStart = min(base["jun"][0],alt["jun"][0])
        #newEnd   = max(a["jun"][-1]+a["len"][-1],b["jun"][-1]+b["len"][-1])
        
        newStart = base["jun"][0]
        newEnd   = base["jun"][-1]+base["len"][-1]
    
        diffL5   = base["jun"][0]-newStart
        diffL3   = newEnd-base["jun"][-1]-base["len"][-1]
        
        newTotL  = base["totalLen"]+diffL5+diffL3
        newPos   = [x + diffL5 for x in base["pos"]]
        newPos[0]= 0                                                  #shift of position starts since the 2nd position. The first position will always be 0.
        
        newStartL= newPos[1]-newPos[0]
        newEndL  = newTotL-newPos[-1]
        
        newJun   = [newStart]+base["jun"][1:]
        newLen   = [newStartL]+base["len"][1:-1]+[newEndL]
        
        newseq5  = base["seq"][base["pos"][0]:base["pos"][1]]
        newseq0  = base["seq"][base["pos"][1]:base["pos"][-1]]
        newseq3  = base["seq"][base["pos"][-1]:]
        aseq5    = a["seq"][a["pos"][0]:a["pos"][1]]
        aseq0    = a["seq"][a["pos"][1]:a["pos"][-1]]
        aseq3    = a["seq"][a["pos"][-1]:]
        bseq5    = b["seq"][b["pos"][0]:b["pos"][1]]
        bseq0    = b["seq"][b["pos"][1]:b["pos"][-1]]
        bseq3    = b["seq"][b["pos"][-1]:]
        
        #if len(aseq5)<len(bseq5) and a["jun"][1]==b["jun"][1]:
        #    newseq5 = bseq5
        #elif len(aseq5)>len(bseq5) and a["jun"][1]==b["jun"][1]:
        #    newseq5 = aseq5
        #if len(aseq3)<len(bseq3) and a["jun"][-1]==b["jun"][-1]:
        #    newseq3 = bseq3
        #elif len(aseq3)>len(bseq3) and a["jun"][-1]==b["jun"][-1]:
        #    newseq3 = aseq3
            
        newseq   = newseq5 + newseq0 + newseq3
        
        newsample = sorted(list(set(a["sample"]+b["sample"])))
        newtissue = sorted(list(set(a["tis"]+b["tis"])))
    
        a["inf"][-3] = ",".join(map(str,newLen))+","                  #block length
        a["inf"][-2] = ",".join(map(str,newPos))+","                  #block relative position
        a["inf"][-1] = ",".join(map(str,newJun))+","                  #block absolute position
        a["inf"][0]  = a["inf"][0]+"|"+b["inf"][0]                    #merging names
        a["inf"][1]  = str(newTotL-int(a["inf"][2]))                  #matches
        a["inf"][3]  = str(max(int(a["inf"][3]),int(b["inf"][3])))    #Number of bases that match but are part of repeats
        a["inf"][4]  = str(max(int(a["inf"][4]),int(b["inf"][4])))    #Number of "N" bases
        a["inf"][5]  = str(max(int(a["inf"][5]),int(b["inf"][5])))    #Number of inserts in query
        a["inf"][6]  = str(max(int(a["inf"][6]),int(b["inf"][6])))    #Number of bases inserted in query
        a["inf"][7]  = str(len(newLen)-1)                             #Number of insert/intron
        a["inf"][8]  = str(max(int(a["inf"][8]),int(b["inf"][8])))    #Number of bases inserted in target
        a["inf"][10] = a["inf"][0]                                    #query name
        a["inf"][11] = str(newTotL)                                   #query size
        a["inf"][13] = str(newTotL)                                   #Alignment end position in query
        a["inf"][16] = map(str,newJun)[0]                             #Alignment start position in target 
        a["inf"][17] = str(newJun[-1]+newLen[-1])                     #Alignment end position in target
        a["inf"][18] = str(len(newLen))                               #Number of block/exon
    
        a["seq"] = newseq
        #if "TTATGGGAGCTATGA" in newseq:
        #    print a["inf"][0],b["inf"][0]
        a["jun"] = newJun
        a["pos"] = newPos
        a["len"] = newLen
        a["totalLen"] = newTotL
        a["sample"] = newsample
        a["tis"] = newtissue
        
        #print diffL5,diffL3
    
    return [a,res]


def ir(a):
    cleanName=[]
    irlist={}                           #e.g.{startPosition:[IR_length1,IR_length2]} 
    for x in a:
        i=0
        while (i+1) < len(a[x]["jun"]):
            left =a[x]["jun"][i]
            right=a[x]["jun"][i+1]
            irlen=right+a[x]["len"][i+1]-left
            if left in irlist:
                if (irlen in irlist[left])==False:
                    irlist[left].append(irlen)
            else:
                irlist[left]=[irlen]
            i+=1
    
    deir={}
    for x in a:
        hit=0
        for y in a[x]["jun"][:-1]:
            for z in irlist[y]:
                if a[x]["len"][a[x]["jun"].index(y)]==z:
                    hit+=1
        if hit==0:
            deir[x]=a[x]
            cleanName.append(x)
                
    return [deir,cleanName]


def merge(v1,v2,transcriptome,a):
    #p=[]
    contig_psl=a+"/pre_TAF1_contigs.psl"
    contig_fa =a+"/pre_TAF1_contigs.fa"
    contig_app=a+"/pre_alltrans_appearance.txt"
    trans_sample=a+"/pre_trans_sample.tab"
    dic={}
    for x in open(contig_psl,"r").readlines():
        t=x[:-1].split("\t")
        name=t[0]
        lng =t[-3]
        jun =t[-1]
        pos =t[-2]
        st  =t[9]
        ch  =t[14]
        dic[name]={}
        dic[name]["jun"]=map(int,jun.split(",")[:-1])
        dic[name]["len"]=map(int,lng.split(",")[:-1])
        dic[name]["pos"]=map(int,pos.split(",")[:-1])
        dic[name]["totalLen"]=int(t[11])
        dic[name]["inf"]=t
        dic[name]["chr"]=ch
        dic[name]["str"]=st
        dic[name]["seq"]=""
        dic[name]["tis"]=""
        dic[name]["sample"]=""
    for x in open(contig_app,"r").readlines():
        t=x[:-1].split("\t")
        name=t[0]
        tissue=sorted(t[1].split(","))
        dic[name]["tis"]=tissue
    for x in open(contig_fa,"r").readlines():
        if x[0]==">":
            name=x[1:-1]
        else:
            dic[name]["seq"]=x[:-1]
    for x in open(trans_sample,"r").readlines():
        t=x[:-1].split("\t")
        name=t[0]
        sample=t[1]
        if name in dic:
            dic[name]["sample"]=sample.split(",")
    
    #print dic
    processed=[]
    for x in dic.keys():
        #if (x in processed)==False:
        for y in dic.keys():       
            #if y!=x and (y in processed) ==False:
            if (y!=x) and (y in dic) and (x in dic):
                test=cmp(dic[x],dic[y])
                dic[x]=test[0]
                if test[1]!=0:
                    #processed.append(y)
                    dic.pop(y,None)
        
        #processed.append(x)
        #p.append(x)
    
    #dic0={}
    #for x in p:
    #    dic0[x]=dic[x]
    
    dic0=dic
    
    deIRres=ir(dic0)
    dic1=deIRres[0]
    q=deIRres[1]
    
    finalID=new(v1,v2,transcriptome,dic1,a)
    
    #print len(dic1),len(finalID)
                    
    f1=open(a+"/TAF1_contigs.psl","w")
    f2=open(a+"/TAF1_contigs.fa","w")
    f3=open(a+"/alltrans_appearance.txt","w")
    f4=open(a+"/trans_sample.tab","w")
    for x in dic1:
        #if x =="NSC_XDP:all:TRINITY_GG_2_c182_g4_i2":
        #    if "TTATGGGAGCTATGA" in dic1[x]["seq"]:
        #        print finalID[x]
        f1.write(finalID[x]+"\t"+"\t".join(dic1[x]["inf"][1:])+"\n")
        f2.write(">"+finalID[x]+"\n")
        f2.write(dic1[x]["seq"]+"\n")
        f3.write(finalID[x]+"\t"+",".join(dic1[x]["tis"])+"\n")
        f4.write(finalID[x]+"\t"+",".join(dic1[x]["sample"])+"\n")
    f1.close()
    f2.close()
    f3.close()
    f4.close()    
    
    getsj(q,finalID,a,"pre_allSamples_SJ.counts.matrix","allSamples_SJ.counts.matrix")
    getsj(q,finalID,a,"pre_allSamples_SJ.fpkm.processed.matrix","allSamples_SJ.fpkm.processed.matrix")
    getsj(q,finalID,a,"pre_allSamples_SJ.normalized.matrix","allSamples_SJ.normalized.matrix")
    
    return 1

merge(70585000,70755500,
      "/data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh37.75_sva_corrected/Homo_sapiens.GRCh37.75.ercc.sva.genePred",
      "/data/talkowski/Samples/XDP/Brain_RNASeq/data/Trinity/RNA_CAP_mergebam/Caudate/all")
    
#merge(70585000,70755500,
      #"/data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh37.75.xdp/Homo_sapiens.GRCh37.75.ercc.sva.genePred",
      #"/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_DN/cpm00")
      
#merge(70585000,70755500,
      #"/data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh37.75.xdp/Homo_sapiens.GRCh37.75.ercc.sva.genePred",
      #"/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_DN_dedup/cpm00")
      
#merge(70585000,70755500,
      #"/data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh37.75.xdp/Homo_sapiens.GRCh37.75.ercc.sva.genePred",
      #"/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_GG/cpm00")
      
#merge(70585000,70755500,
      #"/data/talkowski/tools/ref/RNA-Seq/human/ref_components/GRCh37.75.xdp/Homo_sapiens.GRCh37.75.ercc.sva.genePred",
      #"/data/talkowski/dg520/projects/xdp/output/CellType_up2Jun17_Uniform/merged_transcripts/trinity_GG_dedup/cpm00")
