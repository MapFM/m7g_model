# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 14:01:21 2017

@author: laihongyan
"""
# coding=utf-8
import os

out = open("pse2.txt","w")
f1 = open("Pm7(2).txt", "r+")
for line in f1:
    if line[0] == '>':  # 跳过注释行
        continue
    line = line.upper()  # 将碱基序列转换成大写
    line = line.strip('\n\r')  # 去掉行末换行符
    out.write("1\t")
    for i in range(0, len(line)):
        if line[i] == 'A' or line[i] == 'G':
            x = 1
        else:
            x = 0

        if line[i] == 'A' or line[i] == 'U':
            y = 1
        else:
            y = 0

        if line[i] == 'A' or line[i] == 'C':
            z = 1
        else:
            z = 0
        if i > 1:
            line[0:i-1]
            num = (line[0:i-1].count(line[i],0, i-1))
            d = (num/i)
        else:
            d = 0
        out.write("%d:%d:%d:%f\t" % (x,y,z,d))
    out.write('\n')

f2 = open("Nm7(1).txt", "r+")
for line in f2:
    if line[0] == '>':  # 跳过注释行
        continue
    line = line.upper()  # 将碱基序列转换成大写
    line = line.strip('\n\r')  # 去掉行末换行符
    out.write("2\t")
    for i in range(0, len(line)):
        if line[i] == 'A' or line[i] == 'G':
            x = 1
        else:
            x = 0

        if line[i] == 'A' or line[i] == 'U':
            y = 1
        else:
            y = 0

        if line[i] == 'A' or line[i] == 'C':
            z = 1
        else:
            z = 0
        if i > 1:
            line[0:i-1]
            num = (line[0:i-1].count(line[i],0, i-1))
            d = (num/i)
        else:
            d = 0
        out.write("%d:%d:%d:%f\t" % (x,y,z,d))
    out.write('\n')