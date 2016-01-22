#!/home/andrea/anaconda2/bin/python
import importd as iid

cnlst, cidlst = iid.import_cnames('file3.dat')

# dictionary betweekn kegg and MetaCyc
dic = []
with open('./data/KEGG-Meta.dat','r') as infile:
	for line in infile:
		dic.append(line.strip().split('\t'))

#list of pathway principal ids
pnlst = []
with open('./data/plst_p.dat','r') as infile:
	for line in infile:
		pnlst.append(line.strip().split('\t'))

#list of pathway names
pnnlst = []
with open('./data/plst_n.dat','r') as infile:
	for line in infile:
		pnnlst.append(line.strip())

# create list of IDs in KEGG language
pnlstnew = []
for p in pnlst:
	pnew = []
	for d in dic:
		if d[1] in p and d[1] != '-':
			pnew.append(d[0])	
	pnlstnew.append(pnew)

print pnlstnew
print 'fatto'

# transform the list in program IDs
pidlstnew = []
for p in pnlstnew:
	pnew = []
	for n in p:
		for i,d in enumerate(cnlst):
			if n in d:
				pnew.append(i)
	pidlstnew.append(pnew)
print pidlstnew
print 'fatto'
print pnnlst
summpn = []
summpp = []
for i,p in enumerate(pnnlst):
	if pidlstnew[i] in summpp:
		k = 10000
		for j,x in enumerate(summpp):
			k = j
			if x == pidlstnew[i]:
				break
		summpn[k].append(p)
	else:
		summpn.append([p])
		summpp.append(pidlstnew[i])
print summpp
print summpn
