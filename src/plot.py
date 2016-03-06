#! /usr/bin/env python

import matplotlib.pyplot as plt
import statsmodels.api as sm
import numpy as np
import sys
import re

n = sys.argv[1]
names1 = []
rank1 = []
with open('./results/greedy_l1_rank'+str(n)+'.txt','r') as f1:
	for line in f1:
		l = line.strip().split(',')
		names1.append(l[1])
		rank1.append(abs(float(l[0][1:])))

X = sorted(zip(rank1,names1),reverse=True)
rank1=[x[0] for x in X]
names1=[x[1] for x in X]


names2 = []
rank2 = []
with open('./results/greedy_lda_rank'+str(n)+'.txt','r') as f2:
	for line in f2:
		l = line.strip().split(',')
		names2.append(l[1])
		rank2.append(float(l[0][1:]))

rank3 = [rank2[names2.index(n)] for n in names1[0:20]]

x = rank1[0:20]
y = rank3
print len(x)
print len(y)
m, b = np.polyfit(x, y, 1)
intx = [300000*n for n in x]
y1 = [n-b for n in y]
fig = plt.figure()
ax = fig.add_subplot(111)

ax.scatter(x,y)
ax.plot(x, intx, 'r--')
ax.set_xlim(0,0.004)
ax.set_xlabel('L1-regularized Logistic Regression')
ax.set_ylabel('Linear Discriminant Analysis')
ax.set_title('Prelevance Coefficients of Biomarkers')
for i,z in enumerate(x):
	if x[i] > 0.001 or y[i] > 100:
		if i != 5:
			ax.annotate(re.sub('[^a-zA-Z0-9-]+','',names1[i]),(x[i],y[i]+20.))
		else:
			ax.annotate(re.sub('[^a-zA-Z0-9-]+','',names1[i]),(x[i],y[i]-30.))
#plt.xscale('log')
#plt.yscale('log')
plt.show()
