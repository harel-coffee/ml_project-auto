#!/home/andrea/anaconda2/bin/python
import sys
from sklearn import *
import numpy as np
import cmath as math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sklearn.ensemble import RandomForestClassifier
from analysis import *


def plot_vcurves_irene(Xdata,Ydata,nlst,penlst,param_range):
	"""
	double loo cross validation:
	we leave one sample away and for the remaining set of samples we find the C that maximizes the cross-val score on a loo basis. Then we use that C on the whole set and test it on the point we are interested in. An average over all results for every possible sample left out gives the final result.
	"""
	for pen in penlst:
		score = 1
		counter = 0
		for s1 in cross_validation.LeaveOneOut(len(Ydata)):
			[X1,Y1] = [[Xdata[i] for i in range(len(Xdata)) if i in s1[0]], [Ydata[i] for i in range(len(Ydata)) if i in s1[0]]]
			Cbestscore = 3000
			Cbest = 0
			counter += 1
			print 'sample number'+str(counter)
			for C in param_range:
				Cscore = 1
				for s2 in cross_validation.LeaveOneOut(len(Y1)):
					[X2,Y2] = [[X1[i] for i in range(len(Xdata)) if i in s2[0]], [Y1[i] for i in range(len(Y1)) if i in s2[0]]]
					clf_l1 = linear_model.LogisticRegression(C=C, penalty=pen, tol=1e-3)
					clf_l1.fit(X2, Y2)
					Cscore += int(clf_l1.predict(np.array(X1[s2[1]]).reshape(1,-1))!=[Y1[s2[1]]])
				if Cbestscore > Cscore:
					Cbestscore = Cscore
					Cbest = C
			clf2_l1 = linear_model.LogisticRegression(C=Cbest, penalty=pen, tol=1e-3)
			clf2_l1.fit(X1, Y1)
			score += int(clf2_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
			print int(clf2_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])

		print 'final score: '+str(score)


def plot_vcurves_logistic(Xdata,Ydata,nlst,penlst,param_range):

	"""
	plots validation curves for logistic regression
	"""

	for pen in penlst:

		print 'penalty:\t '+pen
		
		#loo CV
		loo = cross_validation.LeaveOneOut(len(Ydata))
		train_scores1, test_scores1 = learning_curve.validation_curve(linear_model.LogisticRegression(penalty=pen, tol=1e-3),Xdata,Ydata,param_name="C",param_range=param_range,cv=loo)

		train_scores_mean1 = np.mean(train_scores1, axis=1)
		train_scores_std1 = np.std(train_scores1, axis=1)
		test_scores_mean1 = np.mean(test_scores1, axis=1)
		test_scores_std1 = np.std(test_scores1, axis=1)

		print 'LOO -> done'

		#k-fold/bootstrap CV
		bsample = bootstrap(nlst)
		train_scores2, test_scores2 = learning_curve.validation_curve(linear_model.LogisticRegression(penalty=pen, tol=1e-3),Xdata,Ydata,param_name="C",param_range=param_range,cv=10) #substitute with bsample for bootstrap

		train_scores_mean2 = np.mean(train_scores2, axis=1)
		train_scores_std2 = np.std(train_scores2, axis=1)
		test_scores_mean2 = np.mean(test_scores2, axis=1)
		test_scores_std2 = np.std(test_scores2, axis=1)

		print 'k-fold -> done'

		#sparsity calculation
		sparsity_l1_LR = []
		splst = []
		nplst = []
		for C in param_range: 
	
			clf_l1_LR = linear_model.LogisticRegression(C=C, penalty=pen, tol=1e-3)
			clf_l1_LR.fit(Xdata, Ydata)
			coef_l1_LR = clf_l1_LR.coef_.ravel()
			sparsity_l1_LR.append(np.mean(coef_l1_LR == 0) * 100)
			if np.mean(coef_l1_LR != 0)*756 < len(nlst):
				splst.append(C)
			if np.mean(coef_l1_LR == 0) == 1:
				nplst.append(C)

		print 'sparsity -> done'

		#plots
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.semilogx(param_range, train_scores_mean1, label="Training score (loo)", color="r")
		ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
		ax.semilogx(param_range, test_scores_mean1, label="Cross-validation score (loo)",color="g")
		ax.fill_between(param_range, test_scores_mean1 - test_scores_std1,test_scores_mean1 + test_scores_std1, alpha=0.2, color="g")
		ax.semilogx(param_range, train_scores_mean2, label="Training score (bootstrap)", color="k")
		ax.fill_between(param_range, train_scores_mean2 - train_scores_std2, train_scores_mean2 + train_scores_std2, alpha=0.2, color="k")
		ax.semilogx(param_range, test_scores_mean2, label="Cross-validation score (bootstrap)",color="b")
		ax.fill_between(param_range, test_scores_mean2 - test_scores_std2,test_scores_mean2 + test_scores_std2, alpha=0.2, color="b")

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

		ax.set_title("Validation Curve with '"+pen+"' regularized regression and corresponding sparsity")
		fig.savefig("./plots/ValCurveLogit2"+pen+"pen+pca("+str(len(Xdata[1]))+").png")	


def plot_vcurves_SVC(Xdata,Ydata,nlst,penlst,param_range):
	for pen in penlst:

		#loo CV
		print 'penalty:\t '+pen

		loo = cross_validation.LeaveOneOut(len(Ydata))
		train_scores1, test_scores1 = learning_curve.validation_curve(svm.LinearSVC(penalty=pen, dual=False, tol=1e-3),Xdata,Ydata,param_name="C",param_range=param_range,cv=loo)

		train_scores_mean1 = np.mean(train_scores1, axis=1)
		train_scores_std1 = np.std(train_scores1, axis=1)
		test_scores_mean1 = np.mean(test_scores1, axis=1)
		test_scores_std1 = np.std(test_scores1, axis=1)

		print 'LOO -> done'
		
		#kfold CV
		bsample = bootstrap(nlst)
		train_scores2, test_scores2 = learning_curve.validation_curve(svm.LinearSVC(penalty=pen, dual=False, tol=1e-3),Xdata,Ydata,param_name="C",param_range=param_range,cv=10) #substitute with bsample for bootstrap

		train_scores_mean2 = np.mean(train_scores2, axis=1)
		train_scores_std2 = np.std(train_scores2, axis=1)
		test_scores_mean2 = np.mean(test_scores2, axis=1)
		test_scores_std2 = np.std(test_scores2, axis=1)

		print 'k-fold -> done'

		#sparsity calculation
		sparsity_l1_LR = []
		splst = []
		nplst = []
		for C in param_range: 
	
			clf_l1_LR = svm.LinearSVC(C=C, dual=False, penalty=pen, tol=1e-3)
			clf_l1_LR.fit(Xdata, Ydata)
			coef_l1_LR = clf_l1_LR.coef_.ravel()
			sparsity_l1_LR.append(np.mean(coef_l1_LR == 0) * 100)
			if np.mean(coef_l1_LR != 0)*756 < len(nlst):
				splst.append(C)
			if np.mean(coef_l1_LR == 0) == 1:
				nplst.append(C)

		print 'sparsity -> done'

		#plots
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.semilogx(param_range, train_scores_mean1, label="Training score (loo)", color="r")
		ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="r")
		ax.semilogx(param_range, test_scores_mean1, label="Cross-validation score (loo)",color="g")
		ax.fill_between(param_range, test_scores_mean1 - test_scores_std1,test_scores_mean1 + test_scores_std1, alpha=0.2, color="g")
		ax.semilogx(param_range, train_scores_mean2, label="Training score (bootstrap)", color="k")
		ax.fill_between(param_range, train_scores_mean2 - train_scores_std2, train_scores_mean2 + train_scores_std2, alpha=0.2, color="k")
		ax.semilogx(param_range, test_scores_mean2, label="Cross-validation score (bootstrap)",color="b")
		ax.fill_between(param_range, test_scores_mean2 - test_scores_std2,test_scores_mean2 + test_scores_std2, alpha=0.2, color="b")

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


		ax.set_xlabel("$C$")
		ax.set_ylabel("Score")
		ax2.set_ylabel("Sparsity")

		ax.set_title("Validation Curve with '"+pen+"' regularized SVC and corresponding sparsity")
		fig.savefig("./plots/ValCurveSVC"+pen+"pen+pca("+str(len(Xdata[1]))+").png")

def plot_vcurves_RF(Xdata,Ydata,nlst,param_range):
	"""
	plots the cross validation curfe in loo of the random forest classifier
	"""
	loo = cross_validation.LeaveOneOut(len(Ydata))

	test_scores = []
	train_scores = []
	splst = []
	nplst = []

	# analisi anticipata tramite loo per calcolare sparsita (e anche un loo in piu, perche?
	for p in param_range:
		rf = RandomForestClassifier(n_estimators=200, max_depth = p) # works well with 200,5(features)
		rf.fit(Xdata, Ydata)
		#loo cv
		print 'entering the LOO loop for RF'
		looscores = []
		looscores2 = []
		for s in loo:
			X = [Xdata[i] for i in range(len(Xdata)) if i in s[0]]
			Y = [Ydata[i] for i in range(len(Ydata)) if i in s[0]]
			rf.fit(X,Y)
			looscores.append(int(rf.predict(np.array(Xdata[s[1]]).reshape(1,-1))==Ydata[s[1]]))
			looscores2.append(np.mean([int(rf.predict(np.array(Xdata[i]).reshape(1,-1))==Ydata[i]) for i in s[0]]))
		test_scores.append(np.mean(looscores))
		train_scores.append(np.mean(looscores2))
			
	print 'LOO -> done'

	# output description:
	# 
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.semilogx(param_range, train_scores, label="Training score (loo)", color="g")
	ax.semilogx(param_range, test_scores, label="Cross-validation score (loo)",color="b")

	ax.set_xlabel("$p$")
	ax.set_ylabel("Score")

	ax.set_title("Validation Curve for RF")
	fig.savefig("./plots/ValCurveR.png")	

		

def plot_vcurves_pSVC(Xdata,Ydata,Y2data,nlst,penlst,param_range):
	"""
	plots validation curves for support vector machines-based classification in the piecewise-linear case and compares it to the normal one in loo CV
	"""
	for pen in penlst:

		print 'penalty:\t '+pen

		loo = cross_validation.LeaveOneOut(len(Y2data))

		sparsity_l1 = []
		test_scores = []
		splst = []
		nplst = []

		# analisi anticipata tramite loo per calcolare sparsita (e anche un loo in piu, perche?
		for C in param_range: 
			svc_l1 = svm.LinearSVC(C=C, dual=False, penalty=pen, tol=0.01)
			svc_l1.fit(Xdata, Y2data)
			coef_l1 = svc_l1.coef_.ravel()
			sparsity_l1.append(np.mean(coef_l1 == 0) * 100)
			if np.mean(coef_l1 != 0)*756 < len(nlst):
				splst.append(C)
			if np.mean(coef_l1 == 0) == 1:
				nplst.append(C)
			#loo cv
			print 'entering the LOO loop'
			looscores = []
			for s in loo:
				X = [Xdata[i] for i in range(len(Xdata)) if i in s[0]]
				Y = [Y2data[i] for i in range(len(Y2data)) if i in s[0]]
				svc_l1.fit(X,Y)
				looscores.append(mergeclasspred(svc_l1.predict(np.array(Xdata[s[1]]).reshape(1,-1)),Y2data[s[1]]))
			test_scores.append(np.mean(looscores))
				
		print 'LOO-val (standard)'
		train_scores1, test_scores1 = learning_curve.validation_curve(svm.LinearSVC(penalty=pen, dual=False, tol=1e-3),Xdata,Ydata,param_name="C",param_range=param_range,cv=loo)

		train_scores_mean1 = np.mean(train_scores1, axis=1)
		train_scores_std1 = np.std(train_scores1, axis=1)
		test_scores_mean1 = np.mean(test_scores1, axis=1)
		test_scores_std1 = np.std(test_scores1, axis=1)

		print 'LOO -> done'

		# output description:
		# 
		fig = plt.figure(figsize=(10,3))
		ax = fig.add_subplot(111)
		ax.semilogx(param_range, train_scores_mean1, label="Training score (loo)", color="g")
		ax.fill_between(param_range, train_scores_mean1 - train_scores_std1, train_scores_mean1 + train_scores_std1, alpha=0.2, color="g")
		ax.semilogx(param_range, test_scores_mean1, label="Cross-validation score (loo)",color="b")
		ax.fill_between(param_range, test_scores_mean1 - test_scores_std1,test_scores_mean1 + test_scores_std1, alpha=0.2, color="b")
		ax.semilogx(param_range, test_scores, label="Training score (multiloo)", color="r")

		ax2 = ax.twinx()
		ax2.semilogx(param_range,sparsity_l1,color="k")

		try:
			ax.axvline(max(splst),color='r')
		except ValueError:
			ax.axvline(max(param_range),color='r')
		try:
			ax.axvline(max(nplst),color='r')
		except ValueError:
			ax.axvline(min(param_range),color='r')


		ax.set_xlabel("$C$")
		ax.set_ylabel("Score")
		ax2.set_ylabel("Sparsity")

		ax.set_title("Validation Curve with '"+pen+"' regularized SVC and corresponding sparsity")
		fig.savefig("./plots/ValCurvePSVC"+pen+"pen+pca("+str(len(Xdata[1]))+").png")	


def plot_vcurve_loo_l1(Xdata, Ydata, param_range,pen):
	loo = cross_validation.LeaveOneOut(len(Ydata))
	train_scores, test_scores = learning_curve.validation_curve(linear_model.LogisticRegression(penalty=pen, tol=1e-3),Xdata,Ydata,param_name="C",param_range=param_range,cv=10)

	train_scores_mean = np.mean(train_scores, axis=1)
	train_scores_std = np.std(train_scores, axis=1)
	test_scores_mean = np.mean(test_scores, axis=1)
	test_scores_std = np.std(test_scores, axis=1)
	print np.shape(test_scores)
	print sorted([np.mean(t) for t in test_scores])

	plt.title("Validation Curve with regularized regression")
	plt.xlabel("$\lambda$")
	plt.ylabel("Score")
	plt.semilogx(param_range, train_scores_mean, label="Training score", color="r")
	plt.fill_between(param_range, train_scores_mean - train_scores_std, train_scores_mean + train_scores_std, alpha=0.2, color="r")
	plt.semilogx(param_range, test_scores_mean, label="Cross-validation score",color="g")
	plt.fill_between(param_range, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.2, color="g")
	plt.legend(loc="best")
	plt.show()

def plot_vcurve_bootstrap_l1(Xdata, Ydata, nlst):
	bsample = bootstrap(nlst)
	param_range = np.logspace(-7,-6,num=10)
	train_scores, test_scores = learning_curve.validation_curve(linear_model.LogisticRegression(penalty='l1', tol=1e-3),Xdata,Ydata,param_name="C",param_range=param_range,cv=bsample)

	train_scores_mean = np.mean(train_scores, axis=1)
	train_scores_std = np.std(train_scores, axis=1)
	test_scores_mean = np.aray(np.mean(test_scores, axis=1))*0.632 + np.array(np.mean(train_scores, axis=1))*0.368

	test_scores_std = np.std(test_scores, axis=1)
	print np.shape(test_scores)
	print '-----------------------------------------------------\nCI lower bound:\t'+str(sorted([np.mean(t) for t in test_scores])[15])+'\n-----------------------------------------------------'
	
	plt.title("Validation Curve with regularized regression")
	plt.xlabel("$\lambda$")
	plt.ylabel("Score")
	plt.semilogx(param_range, train_scores_mean, label="Training score", color="r")
	plt.fill_between(param_range, train_scores_mean - train_scores_std, train_scores_mean + train_scores_std, alpha=0.2, color="r")
	plt.semilogx(param_range, test_scores_mean, label="Cross-validation score",color="g")
	plt.fill_between(param_range, test_scores_mean - test_scores_std,test_scores_mean + test_scores_std, alpha=0.2, color="g")
	plt.legend(loc="best")
	plt.show()

def cross_vals_score_bootstrap(Xdata, Ydata, nlst, penlst):
	""" this function calculates the empirical lower bound for the confidence interval for the bootstrapped cross validation procedure """
	for pen in penlst:

		scores = []
		for i in range(300):
			print 'progress:\t%.1f%%' % float(i/3.) 
			bsample = bootstrap2(nlst)
			print len(bsample)
			model = linear_model.LogisticRegression(penalty='l1', C=1e-3, tol=1e-3)
			test_scores = cross_validation.cross_val_score(model,Xdata,Ydata,cv=bsample)
			scoreest = np.sum([len(bsample[i][1])*test_scores[i] for i in range(len(test_scores))])/len(Ydata)
			scores.append(scoreest)
		print scores
		sorted(scores)
		print '-----------------------------------------------------\nCI lower bound:\t'+str(sorted(scores)[14])+'\n-----------------------------------------------------'

def bootstrap(nlst):
	""" returns a CV array of bootstrapped training/test sets that together have the same length as the number of predictors"""
	preds = range(len(nlst))
	outlst = []
	while len(outlst) < 300:
		bsample = np.random.choice(preds,size=len(nlst),replace=True)
		outlst.append((np.array([x for x in bsample]),np.array([x for x in range(len(nlst)) if x not in bsample])))
	return outlst

def bootstrap2(nlst):
	""" retorns 1 bootstrap sample of the same size as the training set """
	preds = range(len(nlst))
	bsample = np.random.choice(preds,size=len(nlst),replace=True)
	return [(np.array([x for x in bsample if x != test]),np.array([x for x in bsample if x == test])) for test in preds if test in bsample]

def mergeclasspred(pred,y):
	if (y<4 and pred<4) or (y>3 and pred>3):
		return 1
	else:
		return 0

def plot_vcurves_irene2(Xdata,Ydata,nlst,penlst,param_range,pcadim):
	"""
	double loo cross validation:
	we leave one sample away and for the remaining set of samples we find the C that maximizes the cross-val score on a loo basis. Then we use that C on the whole set and test it on the point we are interested in. An average over all results for every possible sample left out gives the final result.
	"""

	for pen in penlst:	
		cnlst, cidlst = importd.import_cnames('file3.dat')

		with open('./plots/iscores/'+pen+'_'+str(pcadim)+'pca_irene_scores.log','w') as ofile:
			score = 0
			score45 = 0
			counter45 = 0
			score610 = 0
			counter610 = 0
			counter = 0
			#my parameters
			myscore = 0
			myscore45 = 0
			mycounter45 = 0
			myscore610 = 0
			mycounter610 = 0
			mycounter = 0

			myC = 2e-6

			for s1 in cross_validation.LeaveOneOut(len(Ydata)):
				[X1,Y1] = [[Xdata[i] for i in range(len(Xdata)) if i in s1[0]], [Ydata[i] for i in range(len(Ydata)) if i in s1[0]]]
				Cbestscore = 3000
				Cbest = 0
				counter += 1
				name = [nlst[i] for i in range(len(Ydata)) if i in s1[1]][0]

				print 'sample: '+name
				for C in param_range:
					Cscore = 1
					for s2 in cross_validation.LeaveOneOut(len(Y1)):
						[X2,Y2] = [[X1[i] for i in range(len(Xdata)) if i in s2[0]], [Y1[i] for i in range(len(Y1)) if i in s2[0]]]
						clf_l1 = linear_model.LogisticRegression(C=C, penalty=pen, tol=1e-3)
						clf_l1.fit(X2, Y2)
						Cscore += int(clf_l1.predict(np.array(X1[s2[1]]).reshape(1,-1))!=[Y1[s2[1]]])
					if Cbestscore > Cscore:
						Cbestscore = Cscore
						Cbest = C
				#if pcadim == 'all':
					#plot_regpaths2(X1,Y1,name)
					#plot_regpaths3(X1,Y1,name)

				clf2_l1 = linear_model.LogisticRegression(C=Cbest, penalty=pen, tol=1e-3)
				clf2_l1.fit(X1, Y1)
				coef2_l1 = clf2_l1.coef_.ravel()
				sparsity2_l1 = np.mean(coef2_l1 == 0) * 100

				score += int(clf2_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
				print int(clf2_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])

				ofile.write('\n\nscore for '+name+':\t'+str(int(clf2_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]]))+'\nSparsity: '+str(sparsity2_l1)+'\n10 best features:\n'+'\n'.join(map(str,sorted(zip(map(lambda x: round(x, 8), clf2_l1.coef_.ravel()),cidlst),reverse=True))))
				if '10' in name or '6' in name:
					score610 += int(clf2_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
					counter610 += 1
				else:
					score45 += int(clf2_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
					counter45 += 1

				clf3_l1 = linear_model.LogisticRegression(C=myC, penalty=pen, tol=1e-3)
				clf3_l1.fit(X1, Y1)
				myscore += int(clf3_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
				coef3_l1 = clf3_l1.coef_.ravel()
				sparsity3_l1 = np.mean(coef3_l1 == 0) * 100

				print int(clf3_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
				ofile.write('\n\nmyscore for '+name+':\t'+str(int(clf3_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]]))+'\nMysparsity: '+str(sparsity3_l1)+'\n10 best features:\n'+'\n'.join(map(str,sorted(zip(map(lambda x: round(x, 8), clf3_l1.coef_.ravel()),cidlst),reverse=True))))

				if '10' in name or '6' in name:
					myscore610 += int(clf3_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
					mycounter610 += 1
				else:
					myscore45 += int(clf3_l1.predict(np.array(Xdata[s1[1]]).reshape(1,-1))!=[Ydata[s1[1]]])
					mycounter45 += 1


			ofile.write('\n\nFINAL RESULTS:\nScoretot = '+str(score)+'\nCounter = '+str(counter)+'\n\nScore45 = '+str(score45)+'\nCounter45 = '+str(counter45)+'\n\nScore610 = '+str(score610)+'\nCounter610 = '+str(counter610))
				



			print 'final score: '+str(score)
