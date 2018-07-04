#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Conditions:

    OMIM是人类孟德尔遗传疾病在线数据库，记录基因（Gene）条目和表型（Phenotype）条目，每个条目分配唯一的编号，
现在请使用python处理文件omim.txt和mim2gene.txt得到以下统计信息（给出的条目数不是标准答案）
    1).按MIM条目前缀统计（利用文件omim.txt，OMIM号前缀的意义参考)

        MIM Number Prefix	Entry Count
        Gene description *	15530
        Gene and phenotype, combined +	79
        Phenotype description, molecular basis known #	4963
        Phenotype description or locus, molecular basis unknown %	1609
        Other, mainly phenotypes with suspected mendelian basis	1789

    2).按染色体统计基因条目的分布（利用文件omim.txt和mim2gene.txt)基因在X也在Y上面

        chromosome	Entry Count
        chr1	1525
        chr2	969
        chr3	833
        chr4	586
        chr5	696
    
        chrX	719
        chrY	37
        mitochondria	35

"""




import sys
import re
import os

def calNum(infile1, outfile):
    num1,num2,num3,num4,num5=0,0,0,0,0
    dic = {}
    try:
        inF = open(infile1, 'r')
    except IOError:
            sys.exit("open file error")
    outF = open(outfile, 'w')
    lines = inF.readlines()
    for line in lines:
        line = line.strip()
        #print(tmp[0])
        if line.startswith('*'):
            num1+=1
        elif line.startswith('+'):
            num2+=1
        elif line.startswith('#'):
            num3+=1
        elif line.startswith('#'):
            num4+=1
        else :
            num5+=1
    dict = {"Gene description *":str(num1),"Gene and phenotype, combined +":str(num2),
            "Phenotype description, molecular basis known #":str(num3),
            "Phenotype description or locus, molecular basis unknown %":str(num4),
            "Other, mainly phenotypes with suspected mendelian basis":str(num5)}
    for key, value in dict.items():
        outF.write(key+"\t"+value+"\n")
    inF.close()
    outF.close()

def TalFun(inputfile, outputfile):
    arr1 = list(range(1,23))
    arr2 = ["X","Y","mitochondria"]
    arr = arr1+arr2
    num1= 0
    dict = {}
    with open(inputfile, 'r') as inf:
        lines = inf.readlines()
    outf = open(outputfile, 'w')
    for line in lines:
        res = []
        line = line.strip()
        if len(line) > 12:
            tmp = re.split('\s+', line)
            if re.search("(\d+|X|Y|mitochondria)", tmp[2]):
                char = re.search("(\d+|X|Y|mitochondria)", tmp[2])
                res = str(char.group(1))
                for i in range(0,25):
                    a = str(arr[i])
                    if a == res:
                        if re.match('\d+',a):
                            a = "Chr"+a
                            if a not in dict.keys():
                                dict[a] = 0
                            else:
                                dict[a]+=1
                        else:
                            a = a
                            if a not in dict.keys():
                                dict[a] = 0
                            else:
                                dict[a]+=1
                         
            else:
                num1+=1
           # outf.write(tmp[2]+"\t"+str(res)+"\n")
    dict['NA'] = num1
    for key in sorted(dict.keys()):
        outf.write(key+"\t"+str(dict[key])+"\n")
    inf.close()
    outf.close()
def main():
    calNum("omim.txt", "Number.txt")
    TalFun("mim2gene.txt", "Table.txt")
main()
