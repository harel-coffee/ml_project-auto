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
	filename = ''.join(["./data/",organ,"_MS",extension])
	data = pd.read_csv(filename,index_col=0)
	if not isTargeted:
		data = data.transpose()
	#remove rows that are zero
	data = data.loc[:, (data != 0).any(axis=0)]
	print(data)
	nlst = data.index.values[2:]
	Ydata = [1 if 'ctrl' in n else 0 for i,n in enumerate(nlst)]
	Xdata = data.as_matrix()[2:]
	print('Import of "+filename+":\tcomplete')
	return nlst, Xdata, Ydata

def totalimport(organ,isTargeted=True):
	cns, Xdata, Ydata = importdata(organ,isTargeted)
	_, pns = import_cnames(organ,isTargeted)
	return cns, pns, Xdata, Ydata

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


def import_cnames(organ,isTargeted):
	"""
	import the names AND THE KEGG IDS of the chemical elements corresponding to different features of the predictors
	returns
	- a list containing the names of the compounds that have the same mass as the compound ID corresponding to the entry of the array +1
	- a list containing the KEGG-IDs of the compounds that have the same mass as the compound ID corresponding to the entry of the array +1
	"""
	# attention: the number 756 here is correspondent to this experiment!

	if not isTargeted:
		filename = ''.join(["./data/",organ,"_names.csv"])
		data = pd.read_csv(filename,index_col=0)
		KEGGns = data.iloc[:,4]
		cns = data.iloc[:,3]
	else:
		filename = ''.join(["./data/",organ,"_MS_targeted.csv"])
		data = pd.read_csv(filename,index_col=0)
		data = data.loc[:, (data != 0).any(axis=0)]
		KEGGns = data.iloc[:,4]
		cns = data.iloc[:,3]

	return KEGGns, cns

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

def filterd(nlst, Xdata, Ydata, wids=['week_4','week_5','week_6','week_10']):
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
	n, X, Y = importdata(organ)

	print('X:')
	print(X)
	print('-----------------------------------------------------------------------\nY:')
	print(Y)
	print('-----------------------------------------------------------------------\nnlst:')
	print(n)

	print('\n\nData Import for targeted metabolomics\n-----------------------------------------------------------------------\n')
	n, X, Y = importdata('liver',isTargeted=True)
	print('X:')
	print(X)
	print('-----------------------------------------------------------------------\nY:')
	print(Y)
	print('-----------------------------------------------------------------------\nnlst:')
	print(n)

	print('\n--------------------------- Testing Completed Successfully!')

	print('still to do: testing of import_cnames, filterd')

	A, B = import_cnames(organ)

	print(A)
	print(B)
