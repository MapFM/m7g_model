f=open("pse2.txt","r").readlines()
out=open("npf.csv","w")
for line in f:
    a=line.strip("\n").rstrip("\t").split("\t")
    for i in a:
        i = i.split(':')
        for chr in i:
            out.write(chr+',')
    out.write('\n')
