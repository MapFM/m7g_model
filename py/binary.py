f1=open("positive.fasta","r").readlines()
out=open("binary.csv","w")
dic = {
    'A':'1,0,0,0',
    'C':'0,1,0,0',
    'G':'0,0,1,0',
    'U':'0,0,0,1'
}
for i in range(1, len(f1), 2):
    out.write('1'+',')
    for chr in f1[i][:-1]:
        out.write(dic[chr]+',')
    out.write('\n')
f2=open("negative.fasta","r").readlines()
for i in range(1, len(f2), 2):
    print(i)
    out.write('2'+',')
    for chr in f2[i][:-1]:
        out.write(dic[chr]+',')
    out.write('\n')