#coding=utf-8

f=open("1_4_3pse2.txt","r").readlines()


out=open("cn_3_12_csv.csv","w")
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

out.write(','.join(key[0:len(value)]))
out.write('\n')
for m in range(len(all_value)):
	for n in range(len(value)):
		out.write(all_value[m][n]+',')
	out.write('\n')
	
	


