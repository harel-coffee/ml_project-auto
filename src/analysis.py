#! /usr/bin/env python
import sys
from sklearn import *
import numpy as np
import scipy.stats as st
import cmath as math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import importd
import crossval as cv
from scipy import stats

from sklearn.ensemble import RandomForestClassifier
import statsmodels.sandbox.stats.multicomp as mt

def Pvalue(Xdata,Ydata):
	outlst = []
	for j in range(len(Xdata[0])):
		x, pvalue = st.ttest_ind([Xdata[i][j] for i,y in enumerate(Ydata) if y == 1],[Xdata[i][j] for i,y in enumerate(Ydata) if y == 0])
		outlst.append(pvalue)
	return outlst


def qvalues(Xdata,Ydata,compnlst):
	alpha = 0.01
	X1 = np.array([Xdata[i] for i in range(len(Xdata)) if Ydata[i] == 1])
	X0 = np.array([Xdata[i] for i in range(len(Xdata)) if Ydata[i] == 0])
	X1means = [np.mean(X1[:,i]) for i in range(len(Xdata[0]))]
	X0means = [np.mean(X0[:,i]) for i in range(len(Xdata[0]))]

	sigmas = np.array([np.sum([(x[j]-X1means[j])**2 for x in X1])+np.sum([(x[j]-X0means[j])**2 for x in X0]) for j in range(len(Xdata[0]))])/(len(Xdata)-2)
	ses = sigmas*(np.sqrt(1./len(X1)+1./len(X0)))
	print '-----'

	ts = [(X1means[i]-X0means[i])/ses[i] for i in range(len(ses))]
	print ts

	pvalues = [stats.t.sf(np.abs(ts[j]), len(Xdata)-1)*2 for j in range(len(Xdata[0]))]
	print pvalues
	print Pvalue(Xdata,Ydata)
	reject,a,b,c = mt.multipletests(Pvalue(Xdata,Ydata),alpha,method='fdr_bh')
	print reject
	print 'ALPHA = \t'+str(alpha)
	print [c for i,c in enumerate(compnlst) if reject[i]==True]
	

def corr_analysis(Xdata,Ydata, compnlst):
	"""
	correlation and BH analysis fo the features in xdata (whose names are contained in compnlst) and ydata 
	"""
	corrvec = np.array([abs(np.corrcoef(np.array(Xdata)[:,j],Ydata)[0][1]) for j in range(len(Xdata[0]))])
	argsortlst = corrvec.argsort()[::-1]
	#print argsortlst
	s_compnlst = np.array(compnlst)[argsortlst]
	s_corrvec = corrvec[argsortlst]
	return argsortlst, s_compnlst, s_corrvec

def corr_analysis(Xdata,Ydata, compnlst):
	"""
	correlation and BH analysis fo the features in xdata (whose names are contained in compnlst) and ydata 
	"""
	corrvec = np.array([abs(np.corrcoef(np.array(Xdata)[:,j],Ydata)[0][1]) for j in range(len(Xdata[0]))])
	argsortlst = corrvec.argsort()[::-1]
	#print argsortlst
	s_compnlst = np.array(range(len(Xdata[0])))[argsortlst]
	s_corrvec = corrvec[argsortlst]
	return argsortlst, s_compnlst, s_corrvec


def greedy_l1_cv(Xdata,Ydata,addmax,param_range):
	for n in range(1,addmax,10):
		choice = np.array(range(len(Xdata[0])))[argsortlst][:n]
		X = [[x[i] for i in choice] for x in Xdata]
		loo = cross_validation.LeaveOneOut(len(Ydata))
		train_scores1, test_scores1 = learning_curve.validation_curve(linear_model.LogisticRegression(penalty='l1', tol=1e-3),X,Ydata,param_name="C",param_range=param_range,cv=loo)

		train_scores_mean1 = np.mean(train_scores1, axis=1)
		train_scores_std1 = np.std(train_scores1, axis=1)
		test_scores_mean1 = np.mean(test_scores1, axis=1)
		test_scores_std1 = np.std(test_scores1, axis=1)

		sparsity_l1_LR = []
		splst = []
		nplst = []
		for C in param_range: 
	
			clf_l1_LR = linear_model.LogisticRegression(C=C, penalty='l1', tol=1e-3)
			clf_l1_LR.fit(Xdata, Ydata)
			coef_l1_LR = clf_l1_LR.coef_.ravel()
			sparsity_l1_LR.append(np.mean([x == 0 for x in coef_l1_LR]) * 100)
			if np.mean(coef_l1_LR != 0)*n < len(nlst):
				splst.append(C)
			if np.mean(coef_l1_LR == 0) == 1:
				nplst.append(C)


		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.semilogx(param_range, train_scores_mean1, label="Training score (loo)", color="r")
		ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
		ax.semilogx(param_range, test_scores_mean1, label="Cross-validation score (loo)",color="g")
		ax.fill_between(param_range, test_scores_mean1 - test_scores_std1,test_scores_mean1 + test_scores_std1, alpha=0.2, color="g")

		ax2 = ax.twinx()
		ax2.semilogx(param_range,sparsity_l1_LR,color="k")

		try:
			ax.axvline(max(splst),color='r')
		except ValueError:
			ax.axvline(max(param_range),color='r')
		try:
			ax.axvline(max(nplst),color='r')
		except ValueError:
			ax.axvline(min(param_range),color='r')

		ax.legend(loc=0)
		ax2.legend(loc=0)

		ax.set_xlabel("$C$")
		ax.set_ylabel("Score")
		ax2.set_ylabel("Sparsity")

		ax.set_title("Validation Curve with 'l1' regularized regression and corresponding sparsity for n="+str(n)+"features (greedy)")
		fig.savefig("./plots/greedy/ValCurveLogit'l1'pen+greedy("+str(len(X[1]))+").png")

def greedy_l1_boot(Xdata,Ydata,addmax,param_range,cnlst):
	print 'greedy_l1_boot'
	boot = cv.bootstrap(Ydata)
	for n in range(20,addmax,20)+[250,300,350,400,450,500,550,600,700,756]:
		print n
		train_scores = []
		test_scores = []
		sparsities = []
		for b in boot:
		#	print np.shape(np.array(b[0]))
		#	print np.shape(np.array(Xdata))
		#	print np.array(Ydata)
			Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
			 
			argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xboot,Yboot,cnlst)
			choice = np.array(range(len(Xdata[0])))[argsortlst2][:n]
			X2 = [[x[i] for i in choice] for x in Xboot]
			Xtest2 = [[x[i] for i in choice] for x in np.array(Xdata)[b[1]]]
			Ytest = np.array(Ydata)[b[1]]

			sparsity_l1_LR = []
			trscores = []
			tescores = []
			for C in param_range: 
	
				clf_l1_LR = linear_model.LogisticRegression(C=C, penalty='l1', tol=1e-3)
				clf_l1_LR.fit(X2, Yboot)
				trscores.append(clf_l1_LR.score(X2,Yboot))
				tescores.append(clf_l1_LR.score(Xtest2,Ytest))
				
				coef_l1_LR = clf_l1_LR.coef_.ravel()
				sparsity_l1_LR.append(np.mean([x == 0 for x in coef_l1_LR]) * 100)
			train_scores.append(trscores)
			test_scores.append(tescores)
			sparsities.append(sparsity_l1_LR)
		
		train_scores_t = np.array([[train_scores[i][j] for i in range(len(train_scores))] for j in range(len(param_range))])
		test_scores_t = np.array([[test_scores[i][j] for i in range(len(test_scores))] for j in range(len(param_range))])
		train_scores_mean = [np.mean(x) for x in train_scores_t]
		test_scores_mean = np.array([np.mean(x) for x in test_scores_t])*0.632+np.array([np.mean(x) for x in train_scores_t])*0.368
		test_scores_std = [np.sqrt(np.sum(np.array([(test_scores_t[i][j]-test_scores_mean[i])**2 for j in range(len(test_scores[0]))]))/(len(test_scores_t[0])-1)) for i in range(len(test_scores_mean))]


		sparsities_t = np.array([[sparsities[i][j] for i in range(len(sparsities))] for j in range(len(param_range))])	
		sparsities_mean = [np.mean(x) for x in sparsities_t]


		print np.shape(test_scores_mean)
		print np.shape(train_scores_mean)
#		if np.mean(coef_l1_LR != 0)*n < len(nlst):
#			splst.append(C)
#		if np.mean(coef_l1_LR == 0) == 1:
#			nplst.append(C)


		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.semilogx(param_range, train_scores_mean, label="Training score (bootstrap)", color="r")
		#ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
		ax.semilogx(param_range, test_scores_mean, label="Cross-validation score (bootstrap)",color="g")
		ax.fill_between(param_range, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.2, color="g")

		ax2 = ax.twinx()
		ax2.semilogx(param_range,sparsities_mean,color="k")

		#try:
		#	ax.axvline(max(splst),color='r')
		#except ValueError:
		#	ax.axvline(max(param_range),color='r')
		#try:
		#	ax.axvline(max(nplst),color='r')
		#except ValueError:
		#	ax.axvline(min(param_range),color='r')

		ax.legend(loc=0)
		ax2.legend(loc=0)

		ax.set_xlabel("$C$")
		ax.set_ylabel("Score")
		ax2.set_ylabel("Sparsity")

		ax.set_title("Validation Curve with 'l1' regularized regression and corresponding sparsity for n="+str(n)+"features (greedy) in bootstrap")
		fig.savefig("./plots/greedy/ValCurveLogit'l1'pen+greedy("+str(n)+")_bootstrap.png")

def greedy_elnet_boot(Xdata,Ydata,addmax,param_range2,cnlst):
	print 'elastic net'
	Ydata2 = [1 if y==1 else -1 for y in Ydata]
	boot = cv.bootstrap(Ydata2)
	p_range = range(20,addmax,20)
	maxscore = 0
	for n in p_range:
		print n
		train_scores = []
		test_scores = []
		sparsities = []
		for b in boot:
		#	print np.shape(np.array(b[0]))
		#	print np.shape(np.array(Xdata))
		#	print np.array(Ydata)
			Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata2)[b[0]]
			 
			argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xboot,Yboot,cnlst)
			choice = np.array(range(len(Xdata[0])))[argsortlst2][:n]
			X2 = [[x[i] for i in choice] for x in Xboot]
			Xtest2 = [[x[i] for i in choice] for x in np.array(Xdata)[b[1]]]
			Ytest = np.array(Ydata2)[b[1]]

			sparsity_l1_LR = []
			trscores = []
			tescores = []
			for C in param_range2: 
	
				clf_l1_LR = linear_model.ElasticNet(alpha=C, l1_ratio = 0.75, tol=1e-3)
				clf_l1_LR.fit(X2, Yboot)
				trscores.append(clf_l1_LR.score(X2,Yboot))
				tescores.append(clf_l1_LR.score(Xtest2,Ytest))
				
				coef_l1_LR = clf_l1_LR.coef_.ravel()
				sparsity_l1_LR.append(np.mean([x == 0 for x in coef_l1_LR]) * 100)
			train_scores.append(trscores)
			test_scores.append(tescores)
			sparsities.append(sparsity_l1_LR)
		
		print np.shape(np.array(train_scores))
		train_scores_t = np.array([[train_scores[i][j] for i in range(len(train_scores))] for j in range(len(param_range2))])
		test_scores_t = np.array([[test_scores[i][j] for i in range(len(test_scores))] for j in range(len(param_range2))])
		print np.shape(train_scores_t)
		train_scores_mean = [np.mean(x) for x in train_scores_t]
		test_scores_mean = np.array([np.mean(x) for x in test_scores_t])*0.632+np.array([np.mean(x) for x in train_scores_t])*0.368
		test_scores_std = [np.sqrt(np.sum(np.array([(test_scores_t[i][j]-test_scores_mean[i])**2 for j in range(len(test_scores[0]))]))/(len(test_scores_t[0])-1)) for i in range(len(test_scores_mean))]


		sparsities_t = np.array([[sparsities[i][j] for i in range(len(sparsities))] for j in range(len(param_range2))])	
		sparsities_mean = [np.mean(x) for x in sparsities_t]


		print np.shape(test_scores_mean)
		print np.shape(train_scores_mean)
#		if np.mean(coef_l1_LR != 0)*n < len(nlst):
#			splst.append(C)
#		if np.mean(coef_l1_LR == 0) == 1:
#			nplst.append(C)

		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.semilogx(param_range2, train_scores_mean, label="Training score (bootstrap)", color="r")
		#ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
		ax.semilogx(param_range2, test_scores_mean, label="Cross-validation score (bootstrap)",color="g")
		ax.fill_between(param_range2, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.2, color="g")

		ax2 = ax.twinx()
		ax2.semilogx(param_range2,sparsities_mean,color="k")

		#try:
		#	ax.axvline(max(splst),color='r')
		#except ValueError:
		#	ax.axvline(max(param_range),color='r')
		#try:
		#	ax.axvline(max(nplst),color='r')
		#except ValueError:
		#	ax.axvline(min(param_range),color='r')

		ax.legend(loc=0)
		ax2.legend(loc=0)

		ax.set_xlabel("$C$")
		ax.set_ylabel("Score")
		ax2.set_ylabel("Sparsity")

		ax.set_title("Validation Curve with ElasticNet linear regression and corresponding sparsity for n="+str(n)+"features (greedy) in bootstrap")
		fig.savefig("./plots/greedy/ValCurveElNet+greedy("+str(n)+")_bootstrap.png")

		if n == 100:
			True 



def greedy_rf_loo(Xdata,Ydata,addmax,cnlst):
	print 'greedy_rf_loo'

	loo = cross_validation.LeaveOneOut(len(Ydata))
	sparsity_l1_LR = []
	trscores = []
	tescores = []
	tescoresstd = []
	p_range = range(20,addmax,20)+[250,300,350,400,450,500,550,600,700,756]
	sparsities = []
	for n in p_range:
		print n
		train_scores = []
		test_scores = []
		for b in loo:
		#	print np.shape(np.array(b[0]))
		#	print np.shape(np.array(Xdata))
		#	print np.array(Ydata)
			Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
			 
			argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xboot,Yboot,cnlst)
			choice = np.array(range(len(Xdata[0])))[argsortlst2][:n]
			X2 = np.array([[x[i] for i in choice] for x in Xboot])
			Xtest2 = np.array([[x[i] for i in choice] for x in np.array(Xdata)[b[1]]])
			
			Ytest = np.array(Ydata)[b[1]]

			rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
			rf.fit(X2, Yboot)
			train_scores.append(rf.oob_score_)
			test_scores.append(rf.score(Xtest2,Ytest))

		trscores.append(np.mean(train_scores))
		tescores.append(np.mean(test_scores))
		tescoresstd.append(np.std(test_scores))
		sparsities.append((756-n)/756 * 100)
#		if np.mean(coef_l1_LR != 0)*n < len(nlst):
#			splst.append(C)
#		if np.mean(coef_l1_LR == 0) == 1:
#			nplst.append(C)


	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(p_range, trscores, label="Training score (OOB)", color="r")
	#ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
	ax.plot(p_range, tescores, label="Cross-validation score",color="g")
	ax.fill_between(p_range, np.array(tescores) - np.array(tescoresstd),np.array(tescores) + np.array(tescoresstd), alpha=0.2, color="g")

	ax2 = ax.twinx()
	ax2.plot(p_range,sparsities,color="k")

	#try:
	#	ax.axvline(max(splst),color='r')
	#except ValueError:
	#	ax.axvline(max(param_range),color='r')
	#try:
	#	ax.axvline(max(nplst),color='r')
	#except ValueError:
	#	ax.axvline(min(param_range),color='r')

	ax.legend(loc=0)
	ax2.legend(loc=0)

	ax.set_xlabel("$n$")
	ax.set_ylabel("Score")
	ax2.set_ylabel("Sparsity")

	ax.set_title("Validation Curve with greedy random forests and corresponding sparsity for n="+str(n)+"features")
	fig.savefig("./plots/greedy/ValCurveRF+greedy("+str(n)+")_OOB.png")


	train_scores = []
	test_scores = []
	for b in loo:
	#	print np.shape(np.array(b[0]))
	#	print np.shape(np.array(Xdata))
	#	print np.array(Ydata)
		Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
		Xtest = np.array(Xdata)[b[1]]
		Ytest = np.array(Ydata)[b[1]]

		rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
		rf.fit(Xboot, Yboot)
		train_scores.append(rf.score(Xboot,Yboot))
		test_scores.append(rf.score(Xtest,Ytest))

	trscores.append(np.mean(train_scores))
	tescores.append(np.mean(test_scores))
	tescoresstd.append(np.std(test_scores))
	sparsities.append(0)
	#with open('./plots/greedy/greedy_rf_loo.txt','w') as ofile:
	#	p_range = p_range + [756]
	#	ofile.write('\n'.join(['\t'.join([str(p_range[i],str(trscores[i]),str(tescores[i]),str(sparsities[i]))]) for i in range(len(trscores))]))

#		ofile.write('\n\nFeatures of best performance:\n\n')


def greedy_rf_boot(Xdata,Ydata,addmax,cnlst):
	print 'greedy_rf_loo'

	loo = cv.bootstrap(Ydata)
	sparsity_l1_LR = []
	trscores = []
	tescores = []
	tescoresstd = []
	p_range = range(20,addmax,50)+[250,300,400,500,600,700,756]
	sparsities = []
	for n in p_range:
		print n
		train_scores = []
		test_scores = []
		for b in loo:
		#	print np.shape(np.array(b[0]))
		#	print np.shape(np.array(Xdata))
		#	print np.array(Ydata)
			Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]

			argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xboot,Yboot,cnlst)
			choice = np.array(range(len(Xdata[0])))[argsortlst2][:n]
			X2 = np.array([[x[i] for i in choice] for x in Xboot])
			Xtest2 = np.array([[x[i] for i in choice] for x in np.array(Xdata)[b[1]]])
			
			Ytest = np.array(Ydata)[b[1]]

			rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
			rf.fit(X2, Yboot)
			train_scores.append(rf.score(X2,Yboot))
			test_scores.append(rf.score(Xtest2,Ytest))

		trscores.append(np.mean(train_scores))
		tescores.append(np.mean(test_scores))
		tescoresstd.append(np.std(test_scores))
		sparsities.append((756.-n)/756 * 100)
#		if np.mean(coef_l1_LR != 0)*n < len(nlst):
#			splst.append(C)
#		if np.mean(coef_l1_LR == 0) == 1:
#			nplst.append(C)


	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(p_range, trscores, label="Training score (OOB)", color="r")
	#ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
	ax.plot(p_range, tescores, label="Cross-validation score",color="g")
	ax.fill_between(p_range, np.array(tescores) - np.array(tescoresstd),np.array(tescores) + np.array(tescoresstd), alpha=0.2, color="g")

	ax2 = ax.twinx()
	ax2.plot(p_range,sparsities,color="k")

	#try:
	#	ax.axvline(max(splst),color='r')
	#except ValueError:
	#	ax.axvline(max(param_range),color='r')
	#try:
	#	ax.axvline(max(nplst),color='r')
	#except ValueError:
	#	ax.axvline(min(param_range),color='r')

	ax.legend(loc=0)
	ax2.legend(loc=0)

	ax.set_xlabel("$n$")
	ax.set_ylabel("Score")
	ax2.set_ylabel("Sparsity")

	ax.set_title("Validation Curve with greedy random forests and corresponding sparsity for n="+str(n)+"features")
	fig.savefig("./plots/greedy/ValCurveRF+greedy("+str(n)+")_boot(noOOB).png")


	train_scores = []
	test_scores = []
	for b in loo:
	#	print np.shape(np.array(b[0]))
	#	print np.shape(np.array(Xdata))
	#	print np.array(Ydata)
		Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
		Xtest = np.array(Xdata)[b[1]]
		Ytest = np.array(Ydata)[b[1]]

		rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
		rf.fit(Xboot, Yboot)
		train_scores.append(rf.score(Xboot,Yboot))
		test_scores.append(rf.score(Xtest,Ytest))

	trscores.append(np.mean(train_scores))
	tescores.append(np.mean(test_scores))
	tescoresstd.append(np.std(test_scores))
	sparsities.append(0)
	#with open('./plots/greedy/greedy_rf_loo.txt','w') as ofile:
	#	p_range = p_range + [756]
	#	ofile.write('\n'.join(['\t'.join([str(p_range[i],str(trscores[i]),str(tescores[i]),str(sparsities[i]))]) for i in range(len(trscores))]))

#		ofile.write('\n\nFeatures of best performance:\n\n')


def l1_rf_boot(Xdata,Ydata,param_range,cnlst):
	print 'l1_rf_loo'
	loo = cv.bootstrap(Ydata)
	sparsity_l1_LR = []
	trscores = []
	tescores = []
	tescoresstd = []
	sparsities_2 = []
	for n in param_range:
		print n
		train_scores = []
		test_scores = []
		sparsities = []
		for b in loo:
		#	print np.shape(np.array(b[0]))
		#	print np.shape(np.array(Xdata))
		#	print np.array(Ydata)
			Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
			
			clf_l1_LR = linear_model.LogisticRegression(C=n, penalty='l1', tol=1e-3)
			clf_l1_LR.fit(Xboot, Yboot)

			coef_l1_LR = clf_l1_LR.coef_.ravel()
			sparsities.append(np.mean([x == 0 for x in coef_l1_LR]) * 100)

			X2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in Xboot]
			cidlst2 = [cnlst[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]
			 
			Xtest2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in np.array(Xdata)[b[1]]]
			Ytest = np.array(Ydata)[b[1]]

			try:
				rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
				rf.fit(X2, Yboot)
				train_scores.append(rf.score(X2,Yboot))
				test_scores.append(rf.score(Xtest2,Ytest))
			except ValueError:
				train_scores.append(0.5)
				test_scores.append(0.5)


		trscores.append(np.mean(train_scores))
		tescores.append(np.mean(test_scores))
		tescoresstd.append(np.std(test_scores))
		sparsities_2.append(np.mean(sparsities))
#		if np.mean(coef_l1_LR != 0)*n < len(nlst):
#			splst.append(C)
#		if np.mean(coef_l1_LR == 0) == 1:
#			nplst.append(C)


	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.semilogx(param_range, trscores, label="Training score", color="r")
	#ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
	ax.semilogx(param_range, tescores, label="Cross-validation score",color="g")
	ax.fill_between(param_range, np.array(tescores) - np.array(tescoresstd),np.array(tescores) + np.array(tescoresstd), alpha=0.2, color="g")

	ax2 = ax.twinx()
	ax2.semilogx(param_range,sparsities_2,color="k")

	#try:
	#	ax.axvline(max(splst),color='r')
	#except ValueError:
	#	ax.axvline(max(param_range),color='r')
	#try:
	#	ax.axvline(max(nplst),color='r')
	#except ValueError:
	#	ax.axvline(min(param_range),color='r')

	ax.legend(loc=0)

	ax.set_xlabel("$C$")
	ax.set_ylabel("Score")
	ax2.set_ylabel("Sparsity")

	ax.set_title("Validation Curve with l1-logit-preprocessed random forests and corresponding sparsity")
	fig.savefig("./plots/l1/ValCurveRF+l1_boot_2.png")


def l1_rf_l20o(Xdata,Ydata,param_range,cnlst):
	print 'l1_rf_l20'
	loo = l20o(Ydata)
	sparsity_l1_LR = []
	trscores = []
	tescores = []
	tescoresstd = []
	sparsities_2 = []
	splst = []
	nplst = []
	for n in param_range:
		print n
		train_scores = []
		test_scores = []
		sparsities = []
		for b in loo:
		#	print np.shape(np.array(b[0]))
		#	print np.shape(np.array(Xdata))
		#	print np.array(Ydata)
			Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
			
			clf_l1_LR = linear_model.LogisticRegression(C=n, penalty='l1', tol=1e-3)
			clf_l1_LR.fit(Xboot, Yboot)

			coef_l1_LR = clf_l1_LR.coef_.ravel()
			sparsities.append(np.mean([x == 0 for x in coef_l1_LR]) * 100)

			if np.mean([x == 0 for x in coef_l1_LR]) == 0:
				train_scores.append(0.5)
				test_scores.append(0.5)
				continue

			X2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in Xboot]
			cidlst2 = [cnlst[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]
			 
			Xtest2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in np.array(Xdata)[b[1]]]
			Ytest = np.array(Ydata)[b[1]]

			try:
				rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
				rf.fit(X2, Yboot)
				train_scores.append(rf.score(X2,Yboot))
				test_scores.append(rf.score(Xtest2,Ytest))
			except ValueError:
				train_scores.append(0.5)
				test_scores.append(0.5)

		trscores.append(np.mean(train_scores))
		tescores.append(np.mean(test_scores))
		tescoresstd.append(np.std(test_scores))
		sparsities_2.append(np.mean(sparsities))
		
		if np.mean(sparsities)*n < len(nlst):
			splst.append(n)
		if np.mean(sparsities) > 1-1e-3:
			nplst.append(n)


	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.semilogx(param_range, trscores, label="Training score", color="r")
	#ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
	ax.semilogx(param_range, tescores, label="Cross-validation score",color="g")
	ax.fill_between(param_range, np.array(tescores) - np.array(tescoresstd),np.array(tescores) + np.array(tescoresstd), alpha=0.2, color="g")

	ax2 = ax.twinx()
	ax2.semilogx(param_range,sparsities_2,label="Parameter sparsity",color="k")

	try:
		ax.axvline(max(splst),color='r')
	except ValueError:
		ax.axvline(max(param_range),color='r')
	try:
		ax.axvline(max(nplst),color='r')
	except ValueError:
		ax.axvline(min(param_range),color='r')
	ax.set_ylim([0.5,1.1])
	ax.legend(loc=0)

	ax.set_xlabel("$C$")
	ax.set_ylabel("Score")
	ax2.set_ylabel("Sparsity")

	ax.set_title("LOO-CV with l1-logit-preprocessed random forests with sparsity")
	fig.savefig("./plots/l1/ValCurveRF+l1_l20o_3.png")


def l1_rf_loo(Xdata,Ydata,param_range,cnlst):
	print 'l1_rf_loo'
	loo = cross_validation.LeaveOneOut(len(Ydata))
	sparsity_l1_LR = []
	trscores = []
	tescores = []
	tescoresstd = []
	sparsities_2 = []
	splst = []
	nplst = []
	for n in param_range:
		print n
		train_scores = []
		test_scores = []
		sparsities = []
		for b in loo:
		#	print np.shape(np.array(b[0]))
		#	print np.shape(np.array(Xdata))
		#	print np.array(Ydata)
			Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
			
			clf_l1_LR = linear_model.LogisticRegression(C=n, penalty='l1', tol=1e-3)
			clf_l1_LR.fit(Xboot, Yboot)

			coef_l1_LR = clf_l1_LR.coef_.ravel()
			sparsities.append(np.mean([x == 0 for x in coef_l1_LR]) * 100)

			if np.mean([x == 0 for x in coef_l1_LR]) == 0:
				train_scores.append(0.5)
				test_scores.append(0.5)
				continue

			X2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in Xboot]
			cidlst2 = [cnlst[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]
			 
			Xtest2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in np.array(Xdata)[b[1]]]
			Ytest = np.array(Ydata)[b[1]]

			try:
				rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
				rf.fit(X2, Yboot)
				train_scores.append(rf.score(X2,Yboot))
				test_scores.append(rf.score(Xtest2,Ytest))
			except ValueError:
				train_scores.append(0.5)
				test_scores.append(0.5)

		trscores.append(np.mean(train_scores))
		tescores.append(np.mean(test_scores))
		tescoresstd.append(np.std(test_scores))
		sparsities_2.append(np.mean(sparsities))
		
		if np.mean(sparsities)*n < len(nlst):
			splst.append(n)
		if np.mean(sparsities) > 1-1e-3:
			nplst.append(n)


	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.semilogx(param_range, trscores, label="Training score", color="r")
	#ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
	ax.semilogx(param_range, tescores, label="Cross-validation score",color="g")
	ax.fill_between(param_range, np.array(tescores) - np.array(tescoresstd),np.array(tescores) + np.array(tescoresstd), alpha=0.2, color="g")

	ax2 = ax.twinx()
	ax2.semilogx(param_range,sparsities_2,color="k")

	try:
		ax.axvline(max(splst),color='r')
	except ValueError:
		ax.axvline(max(param_range),color='r')
	try:
		ax.axvline(max(nplst),color='r')
	except ValueError:
		ax.axvline(min(param_range),color='r')
	ax.set_ylim([0.5,1.1])
	ax.legend(loc=0)

	ax.set_xlabel("$C$")
	ax.set_ylabel("Score")
	ax2.set_ylabel("Sparsity")

	ax.set_title("LOO-CV with l1-logit-preprocessed random forests with sparsity")
	fig.savefig("./plots/l1/ValCurveRF+l1_loo_2.png")

def l1_LDA_boot(Xdata,Ydata,param_range,cnlst):
	print 'l1_LDA_loo'

	boot = cv.bootstrap(Ydata)
	trscores = []
	tescores = []
	tescoresstd = []
	sparsities_2 = []

	for n in param_range:
		print n
		train_scores = []
		test_scores = []
		sparsities = []
		for b in boot:
		#	print np.shape(np.array(b[0]))
		#	print np.shape(np.array(Xdata))
		#	print np.array(Ydata)
			Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]

			clf_l1_LR = linear_model.LogisticRegression(C=n, penalty='l1', tol=1e-3)
			clf_l1_LR.fit(Xboot, Yboot)

			coef_l1_LR = clf_l1_LR.coef_.ravel()
			print coef_l1_LR, [x == 0 for x in coef_l1_LR]

			sparsities.append(np.mean([x == 0 for x in coef_l1_LR]) * 100)

			X2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in Xboot]
			cidlst2 = [cnlst[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]
			 
			Xtest2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in np.array(Xdata)[b[1]]]
			Ytest = np.array(Ydata)[b[1]]

			clf = discriminant_analysis.LinearDiscriminantAnalysis(solver='eigen',n_components = 2)
			clf.fit(X2,Yboot)

			train_scores.append(clf.score(X2,Yboot))
			test_scores.append(clf.score(Xtest2,Ytest))
	
		trscores.append(np.mean(train_scores))
		tescores.append(np.mean(test_scores))
		tescoresstd.append(np.std(test_scores))
		sparsities_2.append(np.mean(sparsities))

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.semilogx(param_range, trscores, label="Training score", color="r")
	#ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
	ax.semilogx(param_range, tescores, label="Cross-validation score",color="g")
	ax.fill_between(param_range, np.array(tescores) - np.array(tescoresstd),np.array(tescores) + np.array(tescoresstd), alpha=0.2, color="g")

	ax2 = ax.twinx()
	ax2.semilogx(param_range,sparsities_2,color="k")

	#try:
	#	ax.axvline(max(nplst),color='r')
	#except ValueError:
	#	ax.axvline(min(param_range),color='r')

	ax.legend(loc=0)

	ax.set_xlabel("$C$")
	ax.set_ylabel("Score")
	ax2.set_ylabel("Sparsity")

	ax.set_title("Validation Curve with l1-logit-preprocessed LDA and corresponding sparsity")
	fig.savefig("./plots/l1/ValCurveLDA+l1_loo+OOB.png")


def pca(Xdata, nlst, compn):
	"""
	calculates and plots PCA analysis for the data contained in Xdata,
	with names nlst for compn components
	some weeks are selected for the analysis, in order to change these
	weeks the first lines of code must be changed
	"""
	# select some weeks (if week_ is present all weeks are selected)
	w_nlst = [n for i,n in enumerate(nlst) if 'week_10' in nlst[i] or 'week_' in nlst[i]]
	w=[l for i,l in enumerate(Xdata) if 'week_10' in nlst[i] or 'week_' in nlst[i]]
	#separate control and ko
	w_ctrl = [l for i,l in enumerate(w) if 'ctrl' in w_nlst[i]]
	w_ko = [l for i,l in enumerate(w) if 'ko' in w_nlst[i]]

	# define and train the learning algorithm (PCA in this case)
	pca = decomposition.PCA(n_components=compn)
	pca.fit(w)

	Xtrans = pca.transform(w)
	Xtrans_c = pca.transform(w_ctrl)
	Xtrans_k = pca.transform(w_ko)

	# make a plot of the projection on the first three eigendimensions
	try:
		if 'p' in sys.argv[1]:
			xmat = []
			for i in range(3):
				x = list(np.squeeze(Xtrans_k[:,i]))
				x.extend(list(np.squeeze(Xtrans_c[:,i])))
				xmat.append(x)
			xcolorlst = ['r' if i < len(list(np.squeeze(Xtrans_k[:,0]))) else 'b' for i in range(len(x))]
			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')
			ax.scatter(xmat[0],xmat[1],xmat[2],c=xcolorlst)
			fig.saveas('/results/pcatot.png')
	except IndexError:
		print 'no input argument'

	return Xtrans, [0 if 'ko' in n else 1 for n in w_nlst]

def pca2(Xdata, nlst, compn):
	"""
	calculates and plots PCA analysis for the data contained in Xdata,
	with names nlst for compn components
	some weeks are selected for the analysis, in order to change these
	weeks the first lines of code must be changed
	"""
	# select some weeks (if week_ is present all weeks are selected)
	w_nlst = [n for i,n in enumerate(nlst) if 'week_10' in nlst[i] or 'week_' in nlst[i]]
	w=[l for i,l in enumerate(Xdata) if 'week_10' in nlst[i] or 'week_' in nlst[i]]
	#separate control and ko
	w_ctrl = [l for i,l in enumerate(w) if 'ctrl' in w_nlst[i]]
	w_ko = [l for i,l in enumerate(w) if 'ko' in w_nlst[i]]

	# define and train the learning algorithm (PCA in this case)
	pca = decomposition.PCA(n_components=compn)
	pca.fit(w)

	Xtrans = pca.transform(w)
	Xtrans_c = pca.transform(w_ctrl)
	Xtrans_k = pca.transform(w_ko)

	Xtrans2 = pca.inverse_transform(Xtrans)

	# make a plot of the projection on the first three eigendimensions
	try:
		if 'p' in sys.argv[1]:
			xmat = []
			for i in range(3):
				x = list(np.squeeze(Xtrans_k[:,i]))
				x.extend(list(np.squeeze(Xtrans_c[:,i])))
				xmat.append(x)
			xcolorlst = ['r' if i < len(list(np.squeeze(Xtrans_k[:,0]))) else 'b' for i in range(len(x))]
			fig = plt.figure()
			ax = fig.add_subplot(111, projection='3d')
			ax.scatter(xmat[0],xmat[1],xmat[2],c=xcolorlst)
			fig.saveas('/results/pcatot.png')
	except IndexError:
		print 'no input argument'

	return Xtrans2, [0 if 'ko' in n else 1 for n in w_nlst]


#mixed linear model
#minimum t-value

def myFDA(Xdata,Ydata):
	"""
	perform fisher discriminant analysis on the input data
	Xdata is a list of the training data (predictors)
	Ydata is the list of characteristics of each sample
	""" 
	# perform linear discriminant analysis on the given data
	# one can change the solver and the number of components
	print len(Xdata[0])
	clf = discriminant_analysis.LinearDiscriminantAnalysis(solver='eigen',n_components = 2)
	clf.fit(Xdata,Ydata)
	colors = ['red' if i==0 else 'blue' for i in Ydata]

	# transform the data to plot them
	Xplot = np.transpose(clf.transform(Xdata))
	# if we have projected on a single line then transform the data to be able to plot them anyways
	if len(Xplot)==1:
		Xplot = np.append(Xplot,Xplot,axis=0)
	# scatter plot of the transformed data
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.scatter(Xplot[0],Xplot[1],color=colors)
	fig.savefig('./plots/pca('+str(len(Xdata[0]))+')fda.png')
	print '\nScore LDA: '+str(clf.score(Xdata,Ydata))

def myQDA(Xdata,Ydata):
	"""
	prints the score of a QDA on Xdata and Ydata, no parameter needed
	"""
	clf = discriminant_analysis.QuadraticDiscriminantAnalysis()
	clf.fit(Xdata,Ydata)
	print '\nScore QDA: '+str(clf.score(Xdata,Ydata))

def l1_reg(Xdata, Ydata):
	"""
	performs a L1 (lasso) regularized logistic regression on the data that are given as input
	returns (prints) the sparsity of the coefficients of the fitted model and their value
	"""
	# for a range of regularization parameters calculate the sparsity of the fitted l1 regularized model
	nlst, idlst = importd.import_cnames('file3.dat')
	for C in np.logspace(-5,0,num=10): 
	
		clf_l1_LR = linear_model.LogisticRegression(C=C, penalty='l1', tol=1e-3)
		clf_l1_LR.fit(Xdata, Ydata)

		coef_l1_LR = clf_l1_LR.coef_.ravel()
		print len(coef_l1_LR)

		sparsity_l1_LR = np.mean(coef_l1_LR == 0) * 100
		with open('./plots/l1data_C='+str(round(abs(math.log(C,10)),3))+'.dat','w') as outf:
			print("\nC= %.9f" % C)
			print("Sparsity with L1 penalty: %.2f%%" % sparsity_l1_LR)
			print("score with L1 penalty: %.4f" % clf_l1_LR.score(Xdata, Ydata))
			#print sorted(zip(map(lambda x: round(x, 4), clf_l1_LR.coef_.ravel()), range(len(Xdata[0]))), reverse=True)
			outf.write("\nC="+str(C)+'\nSparsity with L1 penalty'+str(sparsity_l1_LR)+'\nscore with L1 penalty: '+str(clf_l1_LR.score(Xdata, Ydata))+'\n\n'+'\n'.join(map(str,sorted(zip(map(lambda x: round(x, 8), clf_l1_LR.coef_.ravel()), idlst), reverse=True))))

def features(idlst, idscores, cnlst):
	""" returns a list containing the scores of the features with nonzero parameter and the respective names """
	return [[idscores[i], cnlst[i]] for i in range(len(idlst))]

def plotfeatures(idlst,idscores,pen,C,Xdata):
	""" plots the regularization coefficients of the different features"""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(range(len(idlst)), idscores, color="r")
	ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train)
	ax.set_xlabel("name$")
	ax.set_ylabel("coefficient")
	ax.set_title("Scores of selected features")
	fig.savefig("./plots/FeaturesScores"+pen+str(C)+"pen+pca("+str(len(Xdata[1]))+").png")	

	

def subtractright(Xdata,nlst,wid1,wid2,C):
	"""
	calculates the set of differences in the data of wid1 and wid2 without mixing the KO and CTRL classes and leaving out one couple of elements each time.
	this is then used to perform a LOO cross-validation of a regularized logistic regression on the sample of the differences of the weeks under consideration
	it returns the average of the scores of the cross-validation, i.e., the value of the validation curve for the penalty C. 
	"""

	w1_nlst = [n for i,n in enumerate(nlst) if wid1 in nlst[i]]
	w1=[l for i,l in enumerate(Xdata) if wid1 in nlst[i]]
	w1_ctrl = [l for i,l in enumerate(w1) if 'ctrl' in w1_nlst[i]]
	w1_ko = [l for i,l in enumerate(w1) if 'ko' in w1_nlst[i]]
	w1_ctrl_nlst = [n for n in w1_nlst if 'ctrl' in n]
	w1_ko_nlst = [n for n in w1_nlst if 'ko' in n]

	w2_nlst = [n for i,n in enumerate(nlst) if wid2 in nlst[i]]
	w2=[l for i,l in enumerate(Xdata) if wid2 in nlst[i]]
	w2_ctrl = [l for i,l in enumerate(w2) if 'ctrl' in w2_nlst[i]]
	w2_ko = [l for i,l in enumerate(w2) if 'ko' in w2_nlst[i]]
	w2_ctrl_nlst = [n for n in w2_nlst if 'ctrl' in n]
	w2_ko_nlst = [n for n in w2_nlst if 'ko' in n]

	scores = []
	for nlo in lolst(nlst,wid1,wid2):
		train_X_ko = np.concatenate([[np.subtract(w2_koi,w1_koj) for i,w2_koi in enumerate(w2_ko) if w2_ko_nlst[i] not in nlo] for j,w1_koj in enumerate(w1_ko) if w1_ko_nlst[j] not in nlo])
		train_X_ctrl = np.concatenate([[np.subtract(w2_ctrli,w1_ctrlj) for i,w2_ctrli in enumerate(w2_ctrl) if w2_ctrl_nlst[i] not in nlo] for j,w1_ctrlj in enumerate(w1_ctrl) if w1_ctrl_nlst[j] not in nlo])
		train_X = np.concatenate((train_X_ko,train_X_ctrl),axis=0)
		train_Y = np.array([0 if i < len(train_X_ko) else 1 for i in range(len(train_X))]) 
		test_X = np.subtract(w2[w2_nlst.index(nlo[1])],w1[w1_nlst.index(nlo[0])])
		test_Y = int('ctrl' in nlo[0])
		model = linear_model.LogisticRegression(penalty='l1', C=C, tol=0.1)
		fit = model.fit(train_X,train_Y)
		scores.append(fit.score([test_X],[test_Y]))
	print '--------------------------------------------------------------------\nRight score for crossing: %.5f' % np.mean(np.array(scores))
	coef_l1_LR = model.coef_.ravel()
	sparsity_l1_LR = np.mean(coef_l1_LR == 0) * 100
	print("\nSparsity with L1 penalty: %.2f%%\n--------------------------------------------------------------------\n\n\n " % sparsity_l1_LR)

	return np.mean(np.array(scores))
		

def lolst(nlst,wid1,wid2):
	"""
	creates a leaveoutlist: a list of all combinations of elements coming from w1 and w2 and from trhe same group (ko and ctrl)
	returns lolst = [(n1k,n2k),...,(n1c,n2c),...,...]
	"""
	nw1c=[n for n in nlst if wid1 in n and 'ctrl' in n]
	nw1k=[n for n in nlst if wid1 in n and 'ko' in n]
	nw2c=[n for n in nlst if wid2 in n and 'ctrl' in n]
	nw2k=[n for n in nlst if wid2 in n and 'ko' in n]
	lolst = []
	for n2 in nw2k:
		lolst.extend([(n1,n2) for n1 in nw1k])
	for n2 in nw2c:
		lolst.extend([(n1,n2) for n1 in nw1c])
	return lolst

def subtract(Xdata,nlst,wid1,wid2):
	"""
	subtraction process, to be commented later
	"""

	w1_nlst = [n for i,n in enumerate(nlst) if wid1 in nlst[i]]
	w1=[l for i,l in enumerate(Xdata) if wid1 in nlst[i]]
	w1_ctrl = [l for i,l in enumerate(w1) if 'ctrl' in w1_nlst[i]]
	w1_ko = [l for i,l in enumerate(w1) if 'ko' in w1_nlst[i]]

	w2_nlst = [n for i,n in enumerate(nlst) if wid2 in nlst[i]]
	w2=[l for i,l in enumerate(Xdata) if wid2 in nlst[i]]
	w2_ctrl = [l for i,l in enumerate(w2) if 'ctrl' in w2_nlst[i]]
	w2_ko = [l for i,l in enumerate(w2) if 'ko' in w2_nlst[i]]

	delta_ctrl = []
	for w1_ctrlj in w1_ctrl:
		delta_ctrl.extend([np.subtract(w2_ctrli,w1_ctrlj) for w2_ctrli in w2_ctrl])

	delta_ko = []
	for w1_koj in w1_ko:
		delta_ko.extend([np.subtract(w2_koi,w1_koj) for w2_koi in w2_ko])

	return delta_ko, delta_ctrl

def plot_regpaths(X,y,param_range):

	eps = 5e-9
	print("Computing regularization path using the lasso...")
	alphas_lasso, coefs_lasso, _ = linear_model.lasso_path(X, y, eps, fit_intercept=True)

	print("Computing regularization path using the positive lasso...")
	alphas_positive_lasso, coefs_positive_lasso, _ = linear_model.lasso_path(X, y, eps, positive=True, fit_intercept=True)

	print("Computing regularization path using the elastic net...")
	alphas_enet, coefs_enet, _ = linear_model.enet_path(X, y, eps=eps, l1_ratio=0.9, fit_intercept=True)

	print("Computing regularization path using the positve elastic net...")
	alphas_positive_enet, coefs_positive_enet, _ = linear_model.enet_path(X, y, eps=eps, l1_ratio=0.9, positive=True, fit_intercept=True)

	#plots
	#ScalarFormatter for ticks

	a = plt.figure(1)
	ax = plt.gca()
	ax.set_color_cycle(2 * ['b', 'r', 'g', 'c', 'k'])
	l1 = plt.plot(-np.log10(alphas_lasso), coefs_lasso.T)
	l2 = plt.plot(-np.log10(alphas_enet), coefs_enet.T, linestyle='--')

	plt.xlabel('-Log(alpha)')
	plt.ylabel('coefficients')
	plt.title('Lasso and Elastic-Net Paths')
	plt.legend((l1[-1], l2[-1]), ('Lasso', 'Elastic-Net'), loc='lower left')
	plt.axis('tight')
	a.savefig('./plots/lasso+elnet/RegPath_Lasso+ElNet.png')

def plot_regpaths2(X,y,name):

	eps = 5e-3
	print("Computing regularization path using the lasso...")
	alphas_lasso, coefs_lasso, _ = linear_model.lasso_path(X, y, eps, fit_intercept=True)

	#print("Computing regularization path using the positive lasso...")
	#alphas_positive_lasso, coefs_positive_lasso, _ = linear_model.lasso_path(X, y, eps, positive=True, fit_intercept=True)
	eps = 5e-10

	print("Computing regularization path using the elastic net...")
	alphas_enet, coefs_enet, _ = linear_model.enet_path(X, y, eps=eps, l1_ratio=0, fit_intercept=True)

	#print("Computing regularization path using the positve elastic net...")
	#alphas_positive_enet, coefs_positive_enet, _ = linear_model.enet_path(X, y, eps=eps, l1_ratio=0.01, positive=True, fit_intercept=True)

	#plots
	#ScalarFormatter for ticks

	a = plt.figure(1,figsize=(9,3))
	a.clf()
	ax = a.add_subplot(111)

	ax.set_color_cycle(2 * ['b', 'r', 'g', 'c', 'k'])

	p2 = ax.plot(-np.log10([x for x in alphas_enet if x < max(alphas_lasso) and x > min(alphas_lasso)]), [coefs_enet.T[i] for i,x in enumerate(alphas_enet) if x < max(alphas_lasso) and x > min(alphas_lasso)], linestyle='--')
	p1 = ax.plot(-np.log10(alphas_lasso), coefs_lasso.T)


	ax.set_xlabel('-Log(alpha)')
	ax.set_ylabel('coefficients')
	ax.set_title('L1 Regularization Paths')
	#ax.legend((p1[-1], p2[-1]), ('L1', 'L2'), loc='lower left')
	#ax.axis('tight')
	a.savefig('./plots/regpaths/'+name.strip().split()[2]+name.strip().split()[0]+'.png')

def plot_regpaths3(X,y,name):

	eps = 5e-1
	print("Computing regularization path using the lasso...")
	alphas_lasso, coefs_lasso, _ = linear_model.lasso_path(X, y, eps, fit_intercept=True)

	#print("Computing regularization path using the positive lasso...")
	#alphas_positive_lasso, coefs_positive_lasso, _ = linear_model.lasso_path(X, y, eps, positive=True, fit_intercept=True)

	print("Computing regularization path using the elastic net...")
	alphas_enet, coefs_enet, _ = linear_model.enet_path(X, y, eps=eps, l1_ratio=0, fit_intercept=True)

	#print("Computing regularization path using the positve elastic net...")
	#alphas_positive_enet, coefs_positive_enet, _ = linear_model.enet_path(X, y, eps=eps, l1_ratio=0.9, positive=True, fit_intercept=True)

	#plots
	#ScalarFormatter for ticks

	a = plt.figure(1,figsize=(9,3))

	ax = plt.gca()
	ax.set_color_cycle(2 * ['b', 'r', 'g', 'c', 'k'])
	#p2 = ax.plot(-np.log10([x for x in alphas_enet if x < max(alphas_lasso) and x > min(alphas_lasso)]), [coefs_enet.T[i] for i,x in enumerate(alphas_enet) if x < max(alphas_lasso) and x > min(alphas_lasso)], linestyle='--')
	p1 = ax.plot(-np.log10(alphas_lasso), coefs_lasso.T)


	ax.set_xlabel('-Log(alpha)')
	ax.set_ylabel('coefficients')
	ax.set_title('L1 and L2 Regularization Paths')
#	ax.legend((p1[-1]), ('L1', 'L2'), loc='lower left')
	ax.axis('tight')
	a.savefig('./plots/regpaths/'+name.strip().split()[2]+name.strip().split()[0]+'_short.png')

def RForest(X,Y,nlst):
	nlst, idlst = importd.import_cnames('file3.dat')
	plst = [2,10,20,50,100,400,None] 
	for p in plst:
		rf = RandomForestClassifier(n_estimators=200, max_depth = p) # works well with 200,5(features)
		rf.fit(X, Y)
		#print "Features sorted by their score:"
		#print sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), range(len(X[0]))), 
        	#     reverse=True)
		with open('./plots/RFdata_P='+str(p)+'.dat','w') as outf:
			#print sorted(zip(map(lambda x: round(x, 4), clf_l1_LR.coef_.ravel()), range(len(Xdata[0]))), reverse=True)
			outf.write('Features extracted by random forest with depth of '+str(p)+'\n'+'\n'.join(map(str,sorted(zip(map(lambda x: round(x, 8), rf.feature_importances_), idlst), reverse=True))))

def RForest_single(X,Y,nlst):
	cnlst, cidlst = importd.import_cnames('file3.dat')
	rf = RandomForestClassifier(n_estimators=300, max_depth = 6,oob_score=True) # works well with 200,5(features)
	rf.fit(X, Y)
	oob = rf.oob_score_
		#print "Features sorted by their score:"
		#print sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), range(len(X[0]))), 
        	#     reverse=True)
	with open('./plots/RFdatatotal.dat','w') as outf:
			#print sorted(zip(map(lambda x: round(x, 4), clf_l1_LR.coef_.ravel()), range(len(Xdata[0]))), reverse=True)
		
		outf.write('Score of the Random Forest algorithm = '+str(oob)+'\n\nFeatures extracted by random forest\n'+'\n'.join(map(str,sorted(zip(map(lambda x: round(x, 8), rf.feature_importances_), cidlst), reverse=True))))

	with open('./plots/RFdataweeks.log','w') as ofile:
		for s1 in cross_validation.LeaveOneOut(len(Y)):
			counter610 = 0
			counter45 = 0
			score610 = 0
			score45 = 0
			name = [nlst[i] for i in range(len(Ydata)) if i in s1[1]]
			[X1,Y1] = [[X[i] for i in range(len(X)) if i in s1[0]], [Y[i] for i in range(len(Y)) if i in s1[0]]]
			rf = RandomForestClassifier(n_estimators=300, max_depth = 6) # works well with 200,5(features)
			rf.fit(X1, Y1)
			print name
			ofile.write('score for '+name[0]+':\t'+str(int(rf.predict(np.array([X[i] for i in range(len(X)) if i in s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]]))+'\nFeatures (top10):\n'+'\n'.join(map(str,sorted(zip(map(lambda x: round(x, 8), rf.feature_importances_), cidlst), reverse=True)[0:10])))

			if '10' in name or '6' in name:
				score610 += int(rf.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
				counter610 += 1
			else:
				score45 += int(rf.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
				counter45 += 1
		ofile.write('\n\nFINAL RESULTS:\nScore45 = '+str(score45)+'\nCounter45 = '+str(counter45)+'\n\nScore610 = '+str(score610)+'\nCounter610 = '+str(counter610))

def RForest_single_greedy(Xdata,Ydata,nlst,m):
	cnlst, cidlst = importd.import_cnames('file3.dat')
	boot = cv.bootstrap(Ydata)
	featimp = np.array([0 for i in range(len(Xdata[0]))])
	train_scores = []
	test_scores = []
	for b in boot:
	#	print np.shape(np.array(b[0]))
	#	print np.shape(np.array(Xdata))
	#	print np.array(Ydata)
		Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
		 
		argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xboot,Yboot,cnlst)
		choice = np.array(range(len(Xdata[0])))[argsortlst2][:m]
		idchoice = [0 for i in range(len(Xdata[0]))]
		for i,c in enumerate(sorted(choice[:])):
			idchoice[c] = i+1
		X2 = np.array([[x[i] for i in choice] for x in Xboot])
		Xtest2 = np.array([[x[i] for i in choice] for x in np.array(Xdata)[b[1]]])
		
		Ytest = np.array(Ydata)[b[1]]

		rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
		rf.fit(X2, Yboot)
		featimp = featimp[:] + np.array([rf.feature_importances_[x-1] if x != 0 else 0 for x in idchoice])
		train_scores.append(rf.score(X2,Yboot))
		test_scores.append(rf.score(Xtest2,Ytest))	#print "Features sorted by their score:"
		#print sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), range(len(X[0]))), 
        	#     reverse=True)
	implst = [f/len(boot) for f in featimp]
	with open('./plots/RFgreedy'+str(m)+'featurelst.txt','w') as ofile:
		ofile.write('\n'.join(map(str,sorted(zip(map(lambda x: round(x, 8), implst), cidlst), reverse=True))))

def l20o(nlst):
	""" returns a CV array of bootstrapped training/test sets that together have the same length as the number of predictors"""
	preds = range(len(nlst))
	outlst = []
	while len(outlst) < 300:
		bsample = np.random.choice(preds,size=92,replace=False)
		outlst.append((np.array([x for x in bsample]),np.array([x for x in range(len(nlst)) if x not in bsample])))
	return outlst

def looo(nlst):
	""" returns a CV array of bootstrapped training/test sets that together have the same length as the number of predictors"""
	preds = range(len(nlst))
	outlst = []
	while len(outlst) < 100:
		bsample = np.random.choice(preds,size=52,replace=False)
		outlst.append((np.array([x for x in bsample]),np.array([x for x in range(len(nlst)) if x not in bsample])))
	return outlst


def RForest_single_l1(Xdata,Ydata,nlst,m):
	cnlst, cidlst = importd.import_cnames('file3.dat')
	boot = l20o(Ydata)
	train_scores = []
	test_scores = []
	tp = 0
	fp = 0
	tn = 0
	fn = 0
	for b in boot:
	#	print np.shape(np.array(b[0]))
	#	print np.shape(np.array(Xdata))
	#	print np.array(Ydata)
		Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
		
		clf_l1_LR = linear_model.LogisticRegression(C=1e-4, penalty='l1', tol=1e-3)
		clf_l1_LR.fit(Xboot, Yboot)

		coef_l1_LR = clf_l1_LR.coef_.ravel()


		if np.mean([x == 0 for x in coef_l1_LR]) == 0:
			train_scores.append(0.5)
			test_scores.append(0.5)
			continue

		X2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in Xboot]
		qidlst2 = [i for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]

		cidlst2 = [cnlst[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]
		 
		Xtest2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in np.array(Xdata)[b[1]]]
		Ytest = np.array(Ydata)[b[1]]

		rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
		rf.fit(X2, Yboot)
		train_scores.append(rf.score(X2,Yboot))
		test_scores.append(rf.score(Xtest2,Ytest))
		predY = rf.predict(Xtest2)
		[(rf.feature_importance_[i],qidlst2[i]) for i in range(len(qidlst2))]
		tp += np.sum([1. if Ytest[i] == 1 and predY[i] == 1 else 0 for i in range(len(predY))])
		fp += np.sum([1. if Ytest[i] == 0 and predY[i] == 1 else 0 for i in range(len(predY))])
		tn += np.sum([1. if Ytest[i] == 0 and predY[i] == 0 else 0 for i in range(len(predY))])
		fn += np.sum([1. if Ytest[i] == 1 and predY[i] == 0 else 0 for i in range(len(predY))])




		#print "Features sorted by their score:"
		#print sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), range(len(X[0]))), 
        	#     reverse=True)
	print tp
	print fp
	print tn
	print fn
#	print 'Rforest_single_l1 loo RESULTS\nAccuracy:\t'+str(np.mean(test_scores))+'\nPrecision:\t'+str(np.mean(precscores))+'\nRecall:\t'+str(np.mean(recallscores))

def RForest_feat_l1(Xdata,Ydata,nlst,m):
	cnlst, cidlst = importd.import_cnames('file3.dat')
	
	clf_l1_LR = linear_model.LogisticRegression(C=1e-3, penalty='l1', tol=1e-3)
	clf_l1_LR.fit(Xdata, Ydata)

	coef_l1_LR = clf_l1_LR.coef_.ravel()

	X2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in Xdata]
	qidlst2 = [i for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]

	cidlst2 = [cnlst[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]
	 

	rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
	rf.fit(X2, Ydata)
	print np.array(cidlst)[qidlst2]
	print '\n'.join(map(str,sorted([(rf.feature_importances_[i],str(cidlst[qidlst2[i]])) for i in range(len(qidlst2))],reverse=True)))

def L1_feat_l1(Xdata,Ydata,nlst,m):
	cnlst, cidlst = importd.import_cnames('file3.dat')
	
	clf_l1_LR = linear_model.LogisticRegression(C=1e-3, penalty='l1', tol=1e-3)
	clf_l1_LR.fit(Xdata, Ydata)

	coef_l1_LR = clf_l1_LR.coef_.ravel()

	X2 = [[x[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0] for x in Xdata]
	qidlst2 = [i for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]

	cidlst2 = [cnlst[i] for i in range(len(Xdata[0])) if coef_l1_LR[i] != 0]
	 

	rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
	rf.fit(X2, Ydata)
	print np.array(cidlst)[qidlst2]
	print '\n'.join(map(str,sorted([(rf.feature_importances_[i],str(cidlst[qidlst2[i]])) for i in range(len(qidlst2))],reverse=True)))



#	print 'Rforest_single_l1 loo RESULTS\nAccuracy:\t'+str(np.mean(test_scores))+'\nPrecision:\t'+str(np.mean(precscores))+'\nRecall:\t'+str(np.mean(recallscores))
def RForest_single_greedy_c(Xdata,Ydata,nlst,m):
	cnlst, cidlst = importd.import_cnames('file3.dat')
	boot = looo(Ydata)
	train_scores = []
	test_scores = []
	print len(list(boot))
	"""
	tp = 0
	fp = 0
	tn = 0
	fn = 0
	for b in boot:
		print tp
	#	print np.shape(np.array(b[0]))
	#	print np.shape(np.array(Xdata))
	#	print np.array(Ydata)
		Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
		

		argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xboot,Yboot,cnlst)
		choice = np.array(range(len(Xdata[0])))[argsortlst2][:80]
		X2 = [[x[i] for i in choice] for x in Xboot]
		Xtest2 = [[x[i] for i in choice] for x in np.array(Xdata)[b[1]]]
		Ytest = np.array(Ydata)[b[1]]

		rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
		rf.fit(X2, Yboot)
		train_scores.append(rf.score(X2,Yboot))
		test_scores.append(rf.score(Xtest2,Ytest))
		predY = rf.predict(Xtest2)
		tp += np.sum([1. if Ytest[i] == 1 and predY[i] == 1 else 0 for i in range(len(predY))])
		fp += np.sum([1. if Ytest[i] == 0 and predY[i] == 1 else 0 for i in range(len(predY))])
		tn += np.sum([1. if Ytest[i] == 0 and predY[i] == 0 else 0 for i in range(len(predY))])
		fn += np.sum([1. if Ytest[i] == 1 and predY[i] == 0 else 0 for i in range(len(predY))])

		#print "Features sorted by their score:"
		#print sorted(zip(map(lambda x: round(x, 4), rf.feature_importances_), range(len(X[0]))), 
        	#     reverse=True)
	print tp
	print fp
	print tn
	print fn

	"""
	rf = RandomForestClassifier(n_estimators=200, max_depth = 6,oob_score=True) # works well with 200,5(features)
	argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xdata,Ydata,cnlst)
	choice = np.array(range(len(Xdata[0])))[argsortlst2][:80]
	X2 = [[x[i] for i in choice] for x in Xdata]

	rf.fit(X2, Ydata)
	featimp = rf.feature_importances_
	importances = sorted([lookfor(featimp,choice,i,cidlst) for i,n in enumerate(cidlst)],reverse=True)
	with open('out_maechler_rf_greedy.txt','w') as fileout:
		fileout.write('\n'.join(map(str,importances)))

#	print 'Rforest_single_l1 loo RESULTS\nAccuracy:\t'+str(np.mean(test_scores))+'\nPrecision:\t'+str(np.mean(precscores))+'\nRecall:\t'+str(np.mean(recallscores))
	return importances



def L1_single_greedy_c(Xdata,Ydata,nlst,m):
	cnlst, cidlst = importd.import_cnames('file3.dat')
	boot = looo(Ydata)
	train_scores = []
	test_scores = []
	print len(list(boot))
	tp = 0
	fp = 0
	tn = 0
	fn = 0
	"""
	for b in boot:
		print tp
		Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
		

		argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xboot,Yboot,cnlst)
		choice = np.array(range(len(Xdata[0])))[argsortlst2][:80]
		X2 = [[x[i] for i in choice] for x in Xboot]
		Xtest2 = [[x[i] for i in choice] for x in np.array(Xdata)[b[1]]]
		Ytest = np.array(Ydata)[b[1]]

		rf = linear_model.LogisticRegression(C=1e-2, penalty='l1', tol=1e-3)
		rf.fit(X2, Yboot)
		train_scores.append(rf.score(X2,Yboot))
		test_scores.append(rf.score(Xtest2,Ytest))
		predY = rf.predict(Xtest2)
		tp += np.sum([1. if Ytest[i] == 1 and predY[i] == 1 else 0 for i in range(len(predY))])
		fp += np.sum([1. if Ytest[i] == 0 and predY[i] == 1 else 0 for i in range(len(predY))])
		tn += np.sum([1. if Ytest[i] == 0 and predY[i] == 0 else 0 for i in range(len(predY))])
		fn += np.sum([1. if Ytest[i] == 1 and predY[i] == 0 else 0 for i in range(len(predY))])
	print tp
	print fp
	print tn
	print fn
	"""

	rf = linear_model.LogisticRegression(C=1e-3, penalty='l1', tol=1e-3)
	argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xdata,Ydata,cnlst)
	choice = np.array(range(len(Xdata[0])))[argsortlst2][:80]
	X2 = [[x[i] for i in choice] for x in Xdata]

	rf.fit(X2, Ydata)
	featimp = map(abs,rf.coef_.ravel())
	print list(s_compnlst)
	print cnlst
	print cnlst[1] in s_compnlst
	importances = sorted([lookfor(featimp,choice,i,cidlst) for i,n in enumerate(cidlst)],reverse=True)
	with open('out_maechler_l1_greedy.txt','w') as fileout:
		fileout.write('\n'.join(map(str,importances)))

	return importances

def lookfor(l1,names, name,cidlst):
	if name in names:
		for i,n in enumerate(names):
			if n == name:
				return (l1[i], cidlst[n])
	else:
		return (0.0, cidlst[name])


def elnet_single_greedy_c(Xdata,Ydata,nlst):
	Ydata = [1 if x==1 else -1 for x in Ydata]
	cnlst, cidlst = importd.import_cnames('file3.dat')
	boot = l20o(Ydata)
	train_scores = []
	test_scores = []
	tp = 0
	fp = 0
	tn = 0
	fn = 0

	for b in boot:
	#	print np.shape(np.array(b[0]))
	#	print np.shape(np.array(Xdata))
	#	print np.array(Ydata)
		Xboot, Yboot = np.array(Xdata)[b[0]],np.array(Ydata)[b[0]]
		 
		argsortlst2, s_compnlst, s_corrvec = corr_analysis(Xboot,Yboot,cnlst)
		choice = np.array(range(len(Xdata[0])))[argsortlst2][:330]
		X2 = [[x[i] for i in choice] for x in Xboot]
		Xtest2 = [[x[i] for i in choice] for x in np.array(Xdata)[b[1]]]
		Ytest = np.array(Ydata)[b[1]]

		clf_l1_LR = linear_model.ElasticNet(alpha=0.4, l1_ratio=0.8, tol=1e-3)
		clf_l1_LR.fit(X2, Yboot)
		predY = clf_l1_LR.predict(Xtest2)
		tp += np.sum([1. if Ytest[i] == 1 and predY[i] == 1 else 0 for i in range(len(predY))])
		fp += np.sum([1. if Ytest[i] == 0 and predY[i] == 1 else 0 for i in range(len(predY))])
		tn += np.sum([1. if Ytest[i] == 0 and predY[i] == 0 else 0 for i in range(len(predY))])
		fn += np.sum([1. if Ytest[i] == 1 and predY[i] == 0 else 0 for i in range(len(predY))])

	print tp
	print fp
	print tn
	print fn





def Foldchange(Xdata,Ydata):
	"""
	returns a list of fold changes in the average of the metabolites at week wid between control and ko mices (ctrl/ko)
	"""
	return [np.mean([Xdata[i][j] for i,y in enumerate(Ydata) if y == 1])/np.mean([Xdata[i][j] for i,y in enumerate(Ydata) if y == 0]) for j in range(len(Xdata[0]))]

def Qvalue(Xdata,Ydata):
	"""
	returns a list of qvalues of the t-test statistics applied on the distributions of each metabolite
	"""
	return True
	# calculate the p stats for each element of ko
	
	# calculate the FWER
	#mt.multipletest

def ExtractNumbers(nlst):
	nlst, Xdata, Ydata = importd.importfile('file.dat')
	nlst, Xdata, Ydata, Y2data = importd.filterd(nlst,Xdata,Ydata,['week_4','week_5','week_6','week_10'])
	a = [0,0,0,0,0,0,0,0]
	for i in Y2data:
		a[i] += 1
	return a


if __name__ == '__main__':
	"""
	
	#define the names of the weeks to be analyzed	
	wids = ['week5','week_6','week_10']
	nlst, Xdata, Ydata = importd.importfile('file.dat')
	cnlst, cidlst = importd.import_cnames('file3.dat')
	# filter the data according to the weeks stored in wids
	nlst, Xdata, Ydata, Y2data = importd.filterd(nlst,Xdata,Ydata,wids)
	#no real need to export pcaY since it is equal to Ydata
	pcaX, pcaY = pca(Xdata, nlst, 300)
	for pca_compn in [10,30,50,70,100,200,300,755]:
		print '------------------------------------------------------------------------\n\n# of PCA components = '+str(pca_compn)+'\n'
		pcaX, pcaY = pca(Xdata, nlst,pca_compn)
		#myFDA(pcaX,pcaY)
	pnlst, pidlst = importd.import_cnames('file3.dat')
	nlst2, Xdata2, Ydata2, Y2data2 = importd.filterd(nlst,Xdata,Ydata,['week_10'])
	foldchanges = Foldchange(Xdata2,Ydata2)
	pvlst = Pvalue(Xdata2,Ydata2)
	print '\n'.join(map(str,[(f, str(pidlst[i])) for i,f in enumerate(foldchanges) if (f > 2 or f < .5) and pvlst[i] < 0.01]))
	print len([(f, str(pidlst[i])) for i,f in enumerate(foldchanges) if (f > 2 or f < .5) and pvlst[i] < 0.01])
	print len([(f, str(pidlst[i])) for i,f in enumerate(foldchanges) if f > 2 or f < .5])
	#print '\n'.join(map(str,sorted([(f, str(pidlst[i])) for i,f in enumerate(foldchanges) if (f > 2 or f < .5) and pvlst[i] < 0.01])))

	#print '\n'.join(map(str,sorted([(f, str(pidlst[i])) for i,f in enumerate(pvlst) if f < 0.01])))

	"""
	
	"""
	for pca_compn in [40,300,755]:
		print '------------------------------------------------------------------------\n\n# of PCA components = '+str(pca_compn)+'\n'
		pcaX, pcaY = pca(Xdata, nlst,pca_compn)
		myFDA(pcaX,pcaY)
		myQDA(pcaX,pcaY)
		l1_reg(pcaX,pcaY)
		cv.plot_vcurve_loo(pcaX,pcaY)
	outlst = cv.bootstrap(nlst)
	cv.plot_vcurve_bootstrap(Xdata,Ydata,nlst)
	
	cv.cross_vals_score_bootstrap(Xdata, Ydata, nlst)
	delta_ko, delta_ctrl = subtract(Xdata, nlst, 'week_4','week_5')
	print len([n for n in nlst if 'ko_week_6' in n])
	print len([n for n in nlst if 'ko_week_10' in n])
	delta_Xdata = np.concatenate((delta_ko,delta_ctrl))
	delta_Ydata = np.array([0 if i < len(delta_ko) else 1 for i in range(len(delta_ko)+len(delta_ctrl))])
	print np.shape(delta_Xdata),np.shape(delta_Ydata)
	l1_reg(delta_Xdata,delta_Ydata)
	cv.plot_vcurve_loo(delta_Xdata,delta_Ydata)
	
	subtractright(Xdata,nlst,'week_5', 'week_6',1e-7)
	"""
	"""
	l1_reg(Xdata, Ydata)
	RForest(Xdata,Ydata,nlst)
	cv.plot_vcurves_RF(Xdata,Ydata,nlst,[2,10,20,50,100,400,None])
	param_range = np.logspace(-8,-3,num=20)
	plot_regpaths(Xdata, Ydata, param_range)
	#cv.plot_vcurves_irene(Xdata,Ydata,nlst,['l1'],param_range)
	cv.plot_vcurves_pSVC(Xdata,Ydata,Y2data,nlst,['l1'],param_range)
	"""



	
	
	penlst = ['l1','l2']
	wids1 = ['week_4','week_5']
	wids2 = ['week_4','week_5','week_6','week_10']
	nlst2, Xdata2, Ydata2 = importd.importfile('file.dat')
	print len(nlst2)
	print len(Ydata2)
	nlst, Xdata, Ydata, Y2data = importd.filterd(nlst2,Xdata2,Ydata2,wids2)
	cnlst, cidlst = importd.import_cnames('file3.dat')
	param_range = np.logspace(-7,-1,num=20)
	print len(Xdata)
	#RForest_single(Xdata,Ydata,nlst)

	#for pca_compn in [40]:
	#	pcaX, pcaY = pca2(Xdata, nlst,pca_compn)
	#	cv.plot_vcurves_irene2(pcaX,pcaY,nlst,penlst,param_range,pca_compn)

	#a = ExtractNumbers(nlst) #outputs the number of samples in each class w4ko, w5ko, ...
	#print a
	#cv.plot_vcurves_irene2(Xdata,Ydata,nlst,penlst,param_range,'all')
		#untransform

	#for pca_compn in [40]:
	#	pcaX, pcaY = pca2(Xdata, nlst,pca_compn)
	#	cv.plot_vcurves_irene2(pcaX,pcaY,nlst,penlst,param_range,pca_compn)
 	#cv.plot_vcurves_logistic(Xdata,Ydata,nlst,penlst,param_range)
	#qvalues(Xdata,Ydata,cidlst)
	#argsortlst, s_compnlst, s_corrvec = corr_analysis(Xdata,Ydata,cidlst)
	#greedy_rf_loo(Xdata,Ydata,200,cidlst)
	#greedy_rf_boot(Xdata,Ydata,200,cidlst)
	#greedy_l1_boot(Xdata,Ydata,200,param_range,cidlst)
#	greedy_elnet_boot(Xdata,Ydata,200,np.logspace(-3,3,num=20),cnlst)
	#l1_rf_l20o(Xdata,Ydata,param_range,cidlst)
	#l1_rf_boot(Xdata,Ydata,param_range,cidlst)
	#l1_LDA_boot(Xdata,Ydata,param_range,cidlst)
	#RForest_single_l1_c(Xdata,Ydata,nlst,1e-3)
	#l1_single_greedy_c(Xdata,Ydata,nlst)
	print '-------------------------------'
	print len(Ydata)
	out1L1 = L1_single_greedy_c(Xdata,Ydata,nlst,1e-3)
	out1Rf = RForest_single_greedy_c(Xdata,Ydata,nlst,1e-3)

	#RForest_single_greedy_c(Xdata,Ydata,nlst,1e-3)
	#elnet_single_greedy_c(Xdata,Ydata,nlst)
	#RForest_feat_l1(Xdata,Ydata,nlst,1)

	nlst, Xdata, Ydata, Y2data = importd.filterd(nlst2,Xdata2,Ydata2,wids1)
	print '-------------------------------'
	print len(Ydata)
	print nlst

	out2L1 = L1_single_greedy_c(Xdata,Ydata,nlst,1e-3)
	out2Rf = RForest_single_greedy_c(Xdata,Ydata,nlst,1e-3)

	print out2Rf

	#plt.plot([i[0] for i in out1L1],[i[0] for i in out2L1],'o')
	#plt.show()
	plt.plot([i[0] for i in out1Rf],[i[0] for i in out2Rf],'o')
	plt.show()






			
