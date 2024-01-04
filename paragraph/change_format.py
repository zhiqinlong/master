# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import vcf 
import re
from Bio import SeqIO
vcf_file=open('D:\工作\工作总结\paragraph\cutesv_sniffles_paragraph.vcf','r')
vcf_reader=vcf.Reader(vcf_file)
ref_file=open('D:\工作\工作总结\paragraph\YS.LG.fasta','r')
aim_vcf=open('D:\工作\工作总结\paragraph\paragraph.txt','w')
ref_dict={}
for record in SeqIO.parse(ref_file,'fasta'):
    ref_dict[record.id]=str(record.seq)
rev=[]
bnd=[]
dup=[]
inv=[]
ins=[]
del1=[]
for record in vcf_reader:
    if record.INFO['SVTYPE']=="TRA" or record.INFO['SVTYPE']=="BND" :
        bnd.append(record.ID)
        continue
    ref=ref_dict[record.CHROM][int(record.POS)-1]
    if record.INFO['END'] > record.POS :
        if record.INFO['SVTYPE']=="DUP" :
            dup.append(record.ID)
            aim_vcf.write(record.CHROM+'\t'+str(record.POS)+'\t'+record.ID+'\t'+ref+'\t'+'<DUP>'+'\t'+'.'+'\t'+'PASS'+'\t'+'SVTYPE='+record.INFO['SVTYPE']+';END='+str(record.INFO['END'])+'\n')
        if record.INFO['SVTYPE']=="INV" :
            inv.append(record.ID)
            aim_vcf.write(record.CHROM+'\t'+str(record.POS)+'\t'+record.ID+'\t'+ref+'\t'+'<INV>'+'\t'+'.'+'\t'+'PASS'+'\t'+'SVTYPE='+record.INFO['SVTYPE']+';END='+str(record.INFO['END'])+'\n')        
        if record.INFO['SVTYPE']=='DEL' :
            del1.append(record.ID)
            aim_vcf.write(record.CHROM+'\t'+str(record.POS)+'\t'+record.ID+'\t'+ref+'\t'+'<DEL>'+'\t'+'.'+'\t'+'PASS'+'\t'+'SVTYPE='+record.INFO['SVTYPE']+';END='+str(record.INFO['END'])+'\n')
    if record.INFO['SVTYPE']== 'INS' :
        ins.append(record.ID)
        seq=ref_dict[record.CHROM][int(record.POS)-1:int(record.POS)+abs(record.INFO['SVLEN'])]
        aim_vcf.write(record.CHROM+'\t'+str(record.POS)+'\t'+record.ID+'\t'+ref+'\t'+seq+'\t'+'.'+'\t'+'PASS'+'\t'+'SVTYPE='+record.INFO['SVTYPE']+'\n')
        
    
len(rev)
len(bnd)
len(inv)
len(ins)
len(del1)              
aim_vcf.close()
ref_file.close()
vcf_file.close()
rev[2]