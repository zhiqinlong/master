#!/usr/bin/python
# -*- coding: utf-8 -*-
import vcf
import re
f1 = open('sniffles_10kb_1bpdis_nonbnd.vcf','r')
f2 = open('sniffles_1bp_seq','r')
output = open('ref_sniffles_10kb_1bpdis_nonbnd.vcf','w')
vcf_reader = vcf.Reader(f1)
ref_seq={}
for eachline in f2:
    eachline=eachline.strip()
    if eachline[0] == '>':
        ids=re.split(':|-',eachline[1:])
        name=ids[0]+':'+str(ids[2]) 
    else:
        ref_seq[str(name)]=eachline
for record in vcf_reader:
    id1=record.CHROM+":"+str(record.POS)
    sv=str(record.ALT[0])
    info="SVTYPE=" + record.INFO['SVTYPE'] + ";END="+str(record.INFO['END'])
    if record.INFO['END'] < record.POS:
        continue
    elif  ">" in  sv or record.INFO['SVTYPE'] == "DUP" or record.INFO['SVTYPE'] == "INV":
        ref=ref_seq[id1]
        record.REF=ref[0]
        if record.INFO['SVTYPE'] == "DUP" :
            output.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(record.CHROM,record.POS,record.ID,record.REF,"<DUP>",".","PASS",info))
        elif record.INFO['SVTYPE'] == "INV":
            output.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(record.CHROM,record.POS,record.ID,record.REF,"<INV>",".","PASS",info))
        else:
            output.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(record.CHROM,record.POS,record.ID,record.REF,sv,".","PASS",info))
    elif record.INFO['SVTYPE'] == "DEL" or  "DEL" in  sv:
        record.REF =str(record.REF)[0]
        output.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(record.CHROM,record.POS,record.ID,record.REF,"<DEL>",".","PASS",info))
    elif record.INFO['SVTYPE'] == "INS":
        alt=str(record.ALT[0])
        record.REF=alt[0]
        output.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n'.format(record.CHROM,record.INFO['END'],record.ID,record.REF,record.ALT[0],".","PASS",info))
f1.close()
f2.close()
output.close()
#####需要uniqstart，因为后面分析发现，有start一样end不一样的变异，这样导致ID有重复CHR:START