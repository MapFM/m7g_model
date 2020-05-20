# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 14:01:21 2017

@author: laihongyan
"""

#coding=utf-8

import os

Twist = [0.063,1.502,0.783,1.071,-1.376,0.063,-1.664,0.783,-0.081,-0.081,0.063,1.502,-1.233,-0.081,-1.376,0.063]
Tilt = [0.502,0.502,0.359,0.215,-1.364,1.077,-1.220,0.359,0.502,0.215,1.077,0.502,-2.368,0.502,-1.364,0.502]
Roll = [0.092,1.195,-0.276,0.827,-1.011,-0.276,-1.379,-0.276,0.092,2.298,-0.276,1.195,-1.379,0.092,-1.011,0.092]
Shift = [1.587,0.126,0.679,-1.019,-0.861,0.560,-0.822,0.679,0.126,-0.348,0.560,0.126,-2.243,0.126,-0.861,1.587]
Slide = [0.111,1.289,-0.241,2.513,-0.623,-0.822,-0.287,-0.241,-0.394,0.646,-0.822,1.289,-1.511,-0.394,-0.623,0.111]
Rise = [-0.109,1.044,-0.623,1.171,-1.254,0.242,-1.389,-0.623,0.711,1.585,0.242,1.044,-1.389,0.711,-1.254,-0.109]
Enthalpy=[1.10,-0.30,0.71,1.44,-0.42,-1.08,0.55,0.71,-1.50,-1.85,-1.08,-0.30,0.51,-0.30,0.71,1.10]
entropy=[0.95,-0.32,0.82,1.42,-0.57,-0.88,0.79,0.82,-1.82,-1.72,-0.88,-0.32,0.27,-0.32,0.82,0.95]
free = [1.56,-0.14,0.07,1.34,0.03,-1.45,-0.29,0.07,-0.28,-1.66,-1.45,-0.14,1.04,-0.28,0.03,1.56]


properties = [Enthalpy, entropy,free]

###将碱基A，C，G，T转换成数字0,1,2,3###
def tran_digital(ch):
	if ch == 'A':
		return 0
	elif ch == 'C':
		return 1
	elif ch == 'G':
		return 2
	else:
		return 3
	
	
###将每次读取的碱基片段转化成十进制标号，如窗口为3时，将AAA到TTT标号为0-63。###
def cal_label(feature):				
	index = 0
	for i in range(0,len(feature)):
		index = 4*index + feature[i]
	return index 

AA=2 #二联体
#以0.05为步长列出从0.05到1
list_w=[]
n=1
while(n>0.95 and n<1.05):
	list_w.append(round(n,1))
	n +=0.1

for para_w in list_w:#以0.05为步长遍历从0.05到1
	for rank in range(4,5):#伪核苷酸的关联从1到50
		for win_size in range(3,4):	#二联体，三联体
		
			filename = str(para_w)+'_'+str(rank)+'_'+str(win_size)
			out = open(filename+"pse2.txt","w")
			
###读取promoter数据,类别标签为1			
			f1 = open("positive.txt","r+")
			for line in f1:	                
				if line[0] == '>':             #跳过注释行
					continue
				line = line.upper()		#将碱基序列转换成大写
				line = line.strip('\n\r')    	#去掉行末换行符  
				number = [0]*len(line)
				
				#将碱基A，C，G，T转换成数字0,1,2,3
				for i in range(0,len(line)):		
					number[i] = tran_digital(line[i])
				
				#求出现频数
				frequence = [0]*(4**win_size)
				for i in range(0,len(number)-win_size+1):   	
					id = cal_label(number[i:i+win_size])
					frequence[id] += 1
				
				#求出现频率
				for i in range(0,len(frequence)):
					frequence[i] /= float(len(number)-win_size+1)
				
				#求Pse特征
				pseall = []
				for m in range(1,rank+1):
					sumpro = 0.0
					for l in range(0,len(properties)):
						for n in range(0,len(number)-AA-m+1):#
							sumpro += properties[l][cal_label(number[n:n+AA])]* properties[l][cal_label(number[n+m:n+m+AA])]
						sumpro /= (len(number)-m-1)
						pseall.append(sumpro)#追加进数组
				sumpse=sum(pseall)
				
					
				###得到并输出所有特征all
				for f in range(0, len(frequence)):
					frequence[f] /=(1+para_w*sumpse)
				for p in range(0, len(pseall)):
					pseall[p] = pseall[p] * para_w / (1+para_w*sumpse)                                                 ###少一个para_w，，，改为: pseall[p] = pseall[p]*para_w / (1+para_w*sumpse) 
				all = frequence + pseall
				out.write("1\t")
				for i in range(len(all)):
					out.write("%d:%.6f\t" % (i+1,all[i]))
				out.write('\n')			#每行结束换行
			f1.close()
			
			
###读取nonpromoter数据，类别标签为2 			
			f2 = open("negative.txt","r+")
			for line in f2:		                
				if line[0] == '>':             #跳过注释行
					continue
				line = line.upper()		#将碱基序列转换成大写
				line = line.strip('\n\r')    	#去掉行末换行符  
				number = [0]*len(line)
				
				#将碱基A，C，G，T转换成数字0,1,2,3
				for i in range(0,len(line)):		
					number[i] = tran_digital(line[i])
				
				#求出现频数
				frequence = [0]*(4**win_size)
				for i in range(0,len(number)-win_size+1):   	
					id = cal_label(number[i:i+win_size])
					frequence[id] += 1
				
				#求出现频率
				for i in range(0,len(frequence)):
					frequence[i] /= float(len(number)-win_size+1)
				
				#求Pse特征
				pseall = []
				for m in range(1,rank+1):
					sumpro = 0.0
					for l in range(0,len(properties)):
						for n in range(0,len(number)-AA-m+1):#
							sumpro += properties[l][cal_label(number[n:n+AA])]* properties[l][cal_label(number[n+m:n+m+AA])]
						sumpro /= (len(number)-m-1)
						pseall.append(sumpro)
				sumpse=sum(pseall)
					
				###得到并输出所有特征all
				for f in range(0, len(frequence)):
					frequence[f] /=(1+para_w*sumpse)
				for p in range(0, len(pseall)):
					pseall[p] = pseall[p] * para_w / (1+para_w*sumpse)
				all = frequence + pseall		
				out.write("2\t")
				for i in range(len(all)):
					out.write("%d:%.6f\t" % (i+293,all[i]))
				out.write('\n')		##每行结束换行		
			f2.close()
			
			out.close()
			#os.system("svm-scale -l 0 -u 1 "+filename+".txt > "+filename+".scale")
			#os.system("gridmy -v 5 "+filename+".scale")
f=open("1_4_3pse2.txt","r").readlines()


out=open("pseKNC.csv","w")
all_value=[]
a=[]
key=[]
key.append('class')

for line in f:
	a=line.strip("\n").split("\t")
	value=[]
	value.append(a[0])
	for i in a[1:-1]:
		value.append(i.split(':')[1])
		key.append(i.split(':')[0])
	all_value.append(value)
	#print(value)

#out.write(','.join(key[0:len(value)]))
#out.write('\n')
for m in range(len(all_value)):
	for n in range(len(value)):
		out.write(all_value[m][n]+',')
	out.write('\n')