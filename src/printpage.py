#!/home/andrea/anaconda2/bin/python

import urllib
import sys

def import_cnames(filename):
	"""
	import the names AND THE KEGG IDS of the chemical elements corresponding to different features of the predictors
	returns
	- cnlst: a list containing the names of the compounds that have the same mass as the compound ID corresponding to the entry of the array +1
	"""
	# attention: the number 756 here is correspondent to this experiment!
	cnlst = []
	cidlst = []
	with open('./data/'+filename,'r') as datafile:
		for i,line in enumerate(datafile):
			l = line.strip().split('\t')
			cnlst.append(l[2].strip().split('; '))
			cidlst.append(l[1].strip().split('; '))
	return cnlst, cidlst

if __name__ == '__main__':
	Knlst, Kidlst = import_cnames('file3.dat')
	Bidlst = []
	with open('./data/KEGG-Meta.dat','w') as ofile:
		for i,M in enumerate(Kidlst):
			for j,Kid in enumerate(M):
				print Kid
				page = urllib.urlopen('http://websvc.biocyc.org/META/foreignid?ids=KEGG:'+Kid)
				page = page.read()
				l = page.strip().split('\t')
				if l[1] == '1':
					Bidlst.append(l[2])
					ofile.write(Kid+'\t'+l[2]+'\n')
				elif l[1] == '0':
					Bidlst.append('-')
					ofile.write(Kid+'\t-\n')

				else:
					print 'CORRESPONDENCE ERROR'
		print 'Correspondence file written'
			
