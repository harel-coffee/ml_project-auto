#! /usr/bin/env python

import sys
import numpy as np
import cmath as math
from functools import reduce
import pandas as pd

# import data
def importdata(organ,isTargeted=False):
	"""
	import the training set from filename
	returns
	- nlst (list of sample names)
	- Xdata (matrix of predictors)
	- Ydata (array of outcomes)
	"""
	if isTargeted:
		extension = '_targeted.csv'
	else:
		extension = '.csv'
	filename = ''.join(["../data/",organ,"_MS",extension])
	data = pd.read_csv(filename,index_col=0)
	if not isTargeted:
		data = data.transpose()
	#remove rows that are zero
	data = data.loc[:, (data != 0).any(axis=0)]
	print(data)
	if not isTargeted:
		pnlst = [i.strip().split(' / ')[-1] for i in data.index.values.tolist()[1:]]
		Ydata = [1 if 'ctrl' in n else 0 for i,n in enumerate(pnlst)]
		Xdata = data.as_matrix()[1:]
		filename2 = ''.join(["../data/",organ,"_names.csv"])
		data2 = pd.read_csv(filename2,index_col=0)
		KEGGns = data2.iloc[:,5]
		cnlst = [i.strip().split('; ')[0] for i in data2.iloc[:,4]]
	else:
		pnlst = data.iloc[:,0].tolist()
		Ydata = [1 if 'ctrl' in n else 0 for i,n in enumerate(pnlst)]
		Xdata = data.as_matrix()[:,2:]
		KEGGns = data.iloc[:,4]
		cnlst = data.columns.values.tolist()[2:]
		pdic = corrnames(organ)
		pnlst = [pdic[i] for i in pnlst]
	print('Import of "+filename+":\tcomplete')
	return pnlst, cnlst, Xdata, Ydata

def corrnames(organ):
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
	if 'week_4' in n or 'week4' in n:
		return 0
	elif 'week_5' in n or 'week5' in n:
		return 1
	elif 'week_6' in n or 'week6' in n:
		return 2
	elif 'week_10' in n or 'week10' in n:
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
	pn, cn, X, Y = importdata('liver',isTargeted=True)
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
	print(len(X))

	corrnames(organ)

	print('\n--------------------------- Import testing Completed Successfully!')
