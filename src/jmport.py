#! /usr/bin/env python

import sys
import numpy as np
import cmath as math
from functools import reduce
import pandas as pd

# import data
def importdata(organ,isTargeted=False,LogConc=None):
	"""
	import the training set from filename
	returns
	- nlst (list of sample names)
	- Xdata (matrix of predictors)
	- Ydata (array of outcomes)
	"""
	#!!!!!!! shitty code ahead
	#remove rows that are zero
	if isTargeted == False:
		extension = '.csv'
		filename = ''.join(["../data/",organ,"_MS",extension])
		data = pd.read_csv(filename,index_col=0)
		data = data.loc[:, (data != 0).any(axis=0)]
		data = data.transpose()
		pnlst = [i.strip().split(' / ')[-1] for i in data.index.values.tolist()[1:]]
		Ydata = [1 if 'ctrl' in n else 0 for i,n in enumerate(pnlst)]
		Xdata = data.as_matrix()[1:]
		filename3 = ''.join(["../data/",organ,"_names.csv"])
		data3 = pd.read_csv(filename3,index_col=0)
		KEGGns = data3.iloc[:,5]
		cnlst = [i.strip().split('; ')[0] for i in data3.iloc[:,4]]
	elif isTargeted == True:
		extension = '_targeted.csv'
		filename = ''.join(["../data/",organ,"_MS",extension])
		data = pd.read_csv(filename,index_col=0)
		data = data.loc[:, (data != 0).any(axis=0)]
		pnlst = data.iloc[:,0].tolist()
		Ydata = [1 if 'ctrl' in n else 0 for i,n in enumerate(pnlst)]
		Xdata = data.as_matrix()[:,2:]
		KEGGns = data.iloc[:,4]
		cnlst = data.columns.values.tolist()[2:]
		pdic = corrnames(organ)
		pnlst = [pdic[i] for i in pnlst]
	else:
		# None

		pnlst1, cnlst1, Xdata1, Ydata1 = importdata(organ,isTargeted=False,LogConc=False)

		pnlst2, cnlst2, Xdata2, Ydata2 = importdata(organ,isTargeted=True,LogConc=False)

		cnlst1.extend(cnlst2)
		X = []
		Ydata = []
		pnlst = []
		for i,p in enumerate(pnlst1):
			for j,q in enumerate(pnlst2):
				if str(q) in p:
					X.append(list(np.concatenate((Xdata1[i],Xdata2[j]),axis=0)))
					Ydata.append(Ydata1[i])
					pnlst.append(p)
					if Ydata1[i] != Ydata2[j]:
						print('FATAL ERROR: wrong identification of\t'+p+'\tand\t'+q)
					break
		Xdata = np.array(X)

	print('Import of "+filename+":\tcomplete')

	if LogConc:
		Xdata = np.log(Xdata)

	elif LogConc == None: #add the logs of the concentratiosn to the Xdata matrix
		pnlst3, cnlst3, Xdata3, Ydata3 = importdata(organ,isTargeted=False,LogConc=True)
		Xdata = np.append(Xdata,Xdata3,1)
		cnlst.extend(['log_'+s for s in cnlst3])

	return pnlst, cnlst, Xdata, Ydata

def corrnames(organ):
	"""
	returns a dictionary correlating the name of the mices in the TM (targeted metabolomics) experiment and the corresponding ID #
	"""
	filename = ''.join(["../data/",organ,"_crossnames.csv"])
	data = pd.read_csv(filename)
	keys = data.iloc[:,2].tolist()
	values = data.iloc[:,1].tolist()
	return dict(zip(keys,values))


# import data
def import_results(filename):
	"""
	Imports the results of analysis.py for statistical analysis and plotting
	"""
	with open(filename,'r') as datafile:

		resultlst = []
		ranklst = []
		scorelst = []
		paramlst = []
		scorestr = ''
		paramstr = ''

		for i,line in enumerate(datafile):
			if line[0] in ['#','\n',' ']:
				continue
			elif line[0] == '-':
				resultlst.append(tmp1)
				ranklst.append(tmp2)
				scorelst.append(scorestr)
				paramlst.append(paramstr)
			elif 'score' in line:
				tmp1 = []
				tmp2 = []
				scorestr = line.strip()
			elif 'parameters' in line:
				paramstr = line.strip()
			else:
				tmp1.append(float(line.strip().split('\t')[0]))
				tmp2.append(int(line.strip().split('\t')[1]))

	return resultlst, ranklst, scorelst, paramlst


def import_cnames(organ,isTargeted=False):
	"""
	import the names AND THE KEGG IDS of the chemical elements corresponding to different features of the predictors
	returns
	- a list containing the names of the compounds that have the same mass as the compound ID corresponding to the entry of the array +1
	- a list containing the KEGG-IDs of the compounds that have the same mass as the compound ID corresponding to the entry of the array +1
	"""
	# attention: the number 756 here is correspondent to this experiment!

	return 0

def casesn(n):
	"""
	returns a number associated to each of the weeks in the range (1; 4)
	"""
	# TODO change to integer switch
	if 'week_4' in n or 'week4' in n or '4w' in n:
		return 0
	elif 'week_5' in n or 'week5' in n or '5w' in n:
		return 1
	elif 'week_6' in n or 'week6' in n or '6w' in n:
		return 2
	elif 'week_10' in n or 'week10' in n or '10w' in n:
		return 3
	else:
		print('error:'+str(n))
		return None

def filterd(nlst, Xdata, Ydata, wids=['week_4','week_5','week_6','week_10','week4','week5','week6','week10']):
	"""
	Filter of the input data according to weeks.
	input:
	- the input data and the names of the weeks you want to select (in an array)
	output:
	- all the data (Xdata, Ydata, nlst) exactly in the same format as before but filtered according to the desired weeks
	- datas of the different families distinguishing between the different weeks: Ydata2 has value from 0 to 7, where the first 4 numbers correspond to "ko" in the corresponding weeks and the second set of 4 elements corresponds to the "ctrl" of the corresponding weeks.
	"""
	y2 = [4*y+casesn(nlst[i]) for i,y in enumerate(Ydata)]
	selectlst = [reduce(lambda x,y: x or y, [wid in n for wid in wids]) for n in nlst]
	return [n for i,n in enumerate(nlst) if selectlst[i]], np.array([x for i,x in enumerate(Xdata) if selectlst[i]]), np.array([y for i,y in enumerate(Ydata) if selectlst[i]]), np.array([y for i,y in enumerate(y2) if selectlst[i]])

def importcombdata(organ,untargeted=True,targeted=False):
	return 0


if __name__ == "__main__":
	organ = 'liver'
	print('Testing...\n\n')
	print('Data Import for normal metabolomics\n-----------------------------------------------------------------------\n')
	pn, cn, X, Y = importdata(organ)

	print('X:')
	print(X)
	print('-----------------------------------------------------------------------\nY:')
	print(Y)
	print('-----------------------------------------------------------------------\npn:')
	print(pn)
	print('-----------------------------------------------------------------------\ncn:')
	print(cn)

	print('\nlength comparison')
	print(len(cn))
	print(len(X[0]))
	print(len(Y))
	print(len(pn))
	print(len(X))

	print('\n\nData Import for targeted metabolomics\n-----------------------------------------------------------------------\n')
	pn, cn, X, Y = importdata('liver',isTargeted=None,LogConc=True)
	print('X:')
	print(X)
	print('-----------------------------------------------------------------------\nY:')
	print(Y)
	print('-----------------------------------------------------------------------\npn:')
	print(pn)
	print('-----------------------------------------------------------------------\ncn:')
	print(cn)

	print('\nlength comparison:')
	print(len(cn))
	print(len(X[0]))
	print(len(Y))
	print(len(pn))
	print(len(X[:,0]))

	corrnames(organ)

	print('\n--------------------------- Import testing Completed Successfully!')
