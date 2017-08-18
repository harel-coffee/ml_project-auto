#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import jmport
from scipy.stats import mannwhitneyu
from operator import itemgetter
from pylab import plot, show, savefig, xlim, figure, \
                hold, ylim, legend, boxplot, setp, axes

#define dictionary for weeks
wdict = {4:0,5:1,6:2,10:3}

# function for setting the colors of the box plots pairs
def setBoxColors(bp):
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    #plt.setp(bp['fliers'][0], markeredgecolor='blue')
    #plt.setp(bp['fliers'][1], markeredgecolor='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    #plt.setp(bp['fliers'][2], markeredgecolor='red')
    #plt.setp(bp['fliers'][3], markeredgecolor='red')
    setp(bp['medians'][1], color='red')


def assignweeks(data,pns):
    out = [[],[],[],[],[],[],[],[]]
    for i,d in enumerate(data):
        out[jmport.casesn(pns[i])+4*int('ko' in pns[i])].append(d)
    return out

def saveplot(data,name,family=1,organ='liver',isTargeted=True):
    """
    function plotting the boxplot data=data[input] of
    metabolite=name[input] (name only used for filename and plot title).
    """
    fig, ax1 = plt.subplots()
    ax1.boxplot(data)
    ax1.set_title('Metabolite: '+name)
    plt.setp(ax1, xticklabels=['ctrl_w4','ctrl_w5','ctrl_w6','ctrl_w10','ko_w4','ko_w5','ko_w6','ko_w10'])
    #plt.setp(xtickNames, rotation=45, fontsize=8)
    fig.savefig('../results/plots/'+name+'_'+organ+'_tar:'+str(isTargeted))

def savegplot(data,name,family=1,organ='liver',isTargeted=True):
    """
    function plotting the groupedboxplot data=data[input] of
    metabolite=name[input] (name only used for filename and plot title)
    where the data is grouped by week and subrouped into ko and ctrl.
    """
    fig, ax = plt.subplots()

    #uncomment the hold command if you want to see the outputted plot
    #hold(True)

    # first boxplot pair (week 4)
    bp = ax.boxplot([data[0],data[4]], positions = [1, 2], widths = 0.6)
    setBoxColors(bp)

    # second boxplot pair (week 5)
    bp = ax.boxplot([data[1],data[5]], positions = [4, 5], widths = 0.6)
    setBoxColors(bp)

    # thrid boxplot pair (week 6)
    bp = ax.boxplot([data[2],data[6]], positions = [7, 8], widths = 0.6)
    setBoxColors(bp)

    # fourth boxplot pair (week10)
    bp = ax.boxplot([data[3],data[7]], positions = [10, 11], widths = 0.6)
    setBoxColors(bp)

    # set axes limits and labels
    xlim(0,12)
    ax.set_xticklabels(['w4', 'w5', 'w6', 'w10'])
    ax.set_xticks([1.5, 4.5, 7.5,10.5])
    ax.set_title('Metabolite: '+name+'   Organ: '+organ)

    # draw temporary red and blue lines and use them to create a legend
    hB, = plot([1,1],'b-')
    hR, = plot([1,1],'r-')
    ax.legend((hB, hR),('Ctrl', 'KO'), loc=2)
    hBs.set_visible(False)
    hR.set_visible(False)

    fig.savefig('../results/plots/'+name+'_'+organ+'_fam'+str(family)+'_tar:'+str(isTargeted)+'.png')
    #fig.show()

def saveps(pns,cns,Xdata,family=1,organ='liver',isTargeted=True):
    """
    function calculating and saving the p values according to mann-whitney
    two-sided U-test for all weeks and all metabolites.
    """
    pdic = [[],[],[],[],[]]
    for i,c in enumerate(cns):
        if 'log' not in c:
            pdic[0].append(c)
    with open('../results/pvals/pvals_fam'+family+'_'+organ+'_tar:'+str(isTargeted)+'.dat','w') as file1:
        for week in range(4):
            print('calculating pvalues for week '+week+'\n')
            wk = wdict[week]
            for i,c in enumerate(cns):
                if 'log' not in c:
                    data = assignweeks(Xdata[:,i],pns)
                    pdic[wk+1].append(mannwhitneyu(data[wk], data[wk+4],alternative="two-sided")[1])
        file1.write('\n'.join([str(i[1])+'\t'+i[0] for i in pdic]))

def savesortedps(week,pns,cns,Xdata,family=1,organ='liver',isTargeted=True):
    """
    function calculating and saving the p values according to mann-whitney
    two-sided U-test for all week=week[input] and all metabolites, sorted
    in decreasing order.
    """
    pdic = []
    with open('../results/pvals/sorted_pvals_week'+str(week)+'fam'+str(family)+'_'+organ+'_tar:'+str(isTargeted)+'.dat','w') as file1:
        for i,c in enumerate(cns):
            if 'log' not in c:
                data = assignweeks(Xdata[:,i],pns)
                print c
                pdic.append((c,mannwhitneyu(data[wdict[week]], data[wdict[week]+4],alternative="two-sided")[1]))
        pdic = sorted(pdic,key=itemgetter(1))
        file1.write('\n'.join([str(i[1])+'\t'+i[0] for i in pdic]))

if __name__ == '__main__':
    organ = 'plasma'
    isTargeted = False
    family=1
    # this need to be changed depending on plasma/liver ---------------------------------------------------------------------------------------------------------------
    if organ == 'liver' and (isTargeted == False or isTargeted == None):
        wids = ['week4','week5','week6','week10']
    elif organ == 'plasma' and isTargeted == False:
        wids = ['week4','week5','week6','week10']
    elif organ == 'liver' and isTargeted == True:
        wids = ['4w','5w','6w','10w']

    pns,cns,Xdata,Ydata = jmport.importdata(organ,family=1,isTargeted=isTargeted)
    _, X, Y, _ = jmport.filterd(pns,Xdata,Ydata,wids)

    print X
    savegplot(assignweeks(Xdata[:,161],pns),'L-rhamnose',family=1,organ=organ) #120 liver and 161 plasma
    savegplot(assignweeks(Xdata[:,162],pns),'L-rhamnose2',family=1,organ=organ) #120 liver and 161 plasma


    for i,c in enumerate(cns):
        print str(i)+'\t'+str(c)
    #     saveplot(assignweeks(Xdata[:,i],pns),c,family=1,organ=organ)

#saveps(pns,cns,Xdata,organ=organ,isTargeted=True)
