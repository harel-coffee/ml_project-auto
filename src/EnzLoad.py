#!/usr/bin/python

import re
from igraph import *
from operator import itemgetter
import itertools
import EnzClass as EC

def decode(organism,dic):
	kdecode(organism,dic)
	print 'kdecode: done'
	cdecode(organism,dic)
	print('cdecode: done\nlen(clst) = '+str(len(dic.clst)))
	qdecode(organism,dic)
	dic.qfillup()
	print('qdecode: done\nlen(qlst) = '+str(len(dic.qlst)))
	rdecode(organism,dic)
	print('rdecode: done\nlen(rlst) = '+str(len(dic.rlst)))
	edecode(organism,dic)
	print('edecode: done\nlen(elst) = '+str(len(dic.elst)))
	print(str(len([e for e in dic.elst if e.activation == True]))+'\tpositive')
	print(str(len([e for e in dic.elst if e.activation == False]))+'\tnegative')
	erdecode(organism,dic)
	print 'regulation check: done'
	dic.create_dgraph()
	threshold = raw_input('threshold (out of dmax='+str(max(dic.dgraph.degree(dic.dgraph.vs[:])))+'): ')
	if threshold in ['','0']:
		filterlst = []
	else:
		filterlst = dic.create_filterlst(int(threshold))
	print 'graphs: done'
	pdecode(organism,dic)
	dic.pflat()
	dic.create_pgraphs(filterlst)
	print('pdecode: done\nlen(plst) = '+str(len(dic.plst)))
	##################################################################################################modify this!!###############################################################3
	colors = [0]
	#mdecode(dic,colors)
	print('mdecode: done')
	print('\nLoading succeeded.')

	print('\nInput/Output reactions:')
	counter = 0
	for r in dic.rlst:
		print r.name,'\t',str(r.slst)
	print('Number of input/output reactions:\t'+str(counter))
	
def kdecode(organism, dic):
	""" decoding of the classes file """
	with open('./organisms/'+organism+'/data/classes.dat','r') as kfile:
		temp = ''
		for line in kfile:
			if line[0] == '#':
				continue
			elif re.search('//',line):
				kstring(temp,dic)
				temp = ''
				continue
			else:
				temp = temp+line+'\n'

def cdecode(organism, dic):
	""" decoding of the compounds file """

	with open('./organisms/'+organism+'/data/compounds.dat','r') as cfile:
		temp = ''
		for line in cfile:
			if line[0] == '#':
				continue
			elif line == '//\n':
				cstring(temp,dic)
				temp = ''
				continue
			else:
				temp = temp+line+'\n'

def qdecode(organism, dic):
	""" decoding of the proteins file """

	with open('./organisms/'+organism+'/data/proteins.dat','r') as qfile:
		temp = ''
		for line in qfile:
			if line[0] == '#':
				continue
			elif line == '//\n':
				qstring(temp,dic)
				temp = ''
				continue
			else:
				temp = temp+line+'\n'


def rdecode(organism, dic):
	""" decoding of the reactions file """

	with open('./organisms/'+organism+'/data/reactions.dat','r') as rfile:
		temp = ''
		for line in rfile:
			if line[0] == '#':
				continue
			elif re.search('//',line):
				rstring(temp,dic)
				temp = ''
				continue
			else:
				temp = temp+line+'\n'


def edecode(organism, dic):
	""" decoding of the enzymes file """

	with open('./organisms/'+organism+'/data/regulation.dat','r') as efile:
		temp = ''
		for line in efile:
			if line[0] == '#':
				continue
			elif re.search('//',line):
				estring(temp,dic)
				temp = ''
				continue
			else:
				temp = temp+line+'\n'

def erdecode(organism, dic):
	""" decoding of the enzrxn file: check of previous steps """

	with open('./organisms/'+organism+'/data/enzrxns.dat','r') as erfile:
		temp = ''
		for line in erfile:
			if line[0] == '#':
				continue
			elif re.search('//',line):
				ercheck(temp,dic)
				temp = ''
				continue
			else:
				temp = temp+line+'\n'


def pdecode(organism, dic):
	""" decoding of the pathways file """

	temp = ''
	allok = False
	firstiteration = True
	while allok == False:
		with open('./organisms/'+organism+'/data/pathways.dat','r') as pfile:
			allok = True
			for line in pfile:
				if line[0:2] == '//':
					#print temp
					fertig = pstring(temp,firstiteration,dic)
					if fertig == False:
						allok = False 
					temp = ''
				else:
					temp = temp+line+'\n'
		firstiteration = False
		print allok
	#self.premdupl()
	
def mdecode(dic,colors):
	""" decoding of the list of inequivalent motifs of n elements given by directg """
	
	with open('./results/motifs/3-motifs_directed.txt','r') as mfile:
		for line in mfile:
			mstring(line,colors,dic)

def kstring(string,dic):
	""" splits and stores the strings coming from the classes file """

	kname = ''
	cmname = []
	types = []
	
	lines = string.strip().split('\n')
	for i in lines:
		line = i.strip().split(' - ')
		if line[0] == 'UNIQUE-ID':
			kname = line[1]

		if line[0] == 'COMMON-NAME':
			cmname.append(line[1])

		if line[0] == 'TYPES':
			types.append(line[1])

	k = dic.ksearch(kname)
	if k == None:
		k = EC.Class(kname)
		dic.appendk(k)

	for tname in types:
		t = dic.ksearch(tname)
		if t == None:
			t = EC.Class(tname)
			dic.appendk(t)
		t.addchild(k)
		k.addfather(t)

	k.setcommname(cmname)


def cstring(string,dic):
	""" splits and stores the strings coming from the compounds file """

	name = ''
	commname = []
	types = []
	reglst = []
	smiles = ''
	weight = 0
	
	lines = string.strip().split('\n')
	for i in lines:
		line = i.strip().split(' - ')
		if line[0] == 'UNIQUE-ID':
			name = line[1]

		if line[0] == 'COMMON-NAME':
			commname.append(line[1])

		if line[0] == 'TYPES':
			k = dic.ksearch(line[1])
			if k == None:
				k = EC.Class(line[1])
				dic.appendk(k)
			types.append(k)

		if line[0] == 'SMILES':
			smiles = line[1]

		if line[0] == 'WEIGHT':
			weight = int(line[1])

		if line[0] == 'REGULATES':
			reglst.append(line[1])

	if name == '' or types == [] or smiles == '':
		print('- missing information for compound '+name+'\n')

	if dic.csearch(name) == None:
		c = EC.Compound(name, commname, types, smiles, weight, reglst)
		dic.appendc(c)
	else: 
		#error: compound already present
		print('- compound '+name+' already existing in database at position'+str(dic.csearch(name))+'\nThe input is:\n'+string)
	for k in c.typ:
		k.addg(c)


def qstring(string,dic):
	""" splits and stores the strings coming from the proteins file """

	name = ''
	commname = []
	types = []
	comp = []
	compof = []
	reg = []
	enz = []
	
	lines = string.strip().split('\n')
	for i in lines:
		line = i.strip().split(' - ')
		if line[0] == 'UNIQUE-ID':
			name = line[1]

		elif line[0] == 'COMMON-NAME':
			commname.append(line[1])

		elif line[0] == 'TYPES':
			k = dic.ksearch(line[1])
			if k == None:
				k = EC.Class(line[1])
				dic.appendk(k)
			types.append(k)

		elif line[0] == 'COMPONENTS':
			comp.append(line[1])

		elif line[0] == '^COEFFICIENT':
			for i in range(int(line[1])-1):
				comp.append(comp[len(comp)-1])

		elif line[0] == 'REGULATES':
			reg.append(line[1])

		elif line[0] == 'COMPONENT-OF':
			compof.append(line[1])

		elif line[0] == 'CATALYZES':
			enz.append(line[1])


	if name == '' or types == []:
		print('- missing information for compound '+name+'\n')

	if dic.qsearch(name) == None:
		q = EC.Protein(name, commname, types, comp, compof, reg, enz)
		dic.appendq(q)
	else: 
		#error: compound already present
		print('- protein '+name+' already existing in database at position'+str(dic.qsearch(name))+'\n')

	for k in q.typ:
		k.addg(q)


def rstring(string,dic):
	""" splits and stores the strings coming from the reactions file """

	rname = ''
	typ = []
	commname = []
	enames = []
	innames = []
	outnames = []
	slst = []
	tlst = []
	nobalance = False
	inpath = []
	reversible = True
	lines = string.strip().split('\n')
	reverse = False
	for i in lines:
		line = i.strip().split(' - ')
		if line[0] == 'UNIQUE-ID':
			rname = line[1]
		if line[0] == 'REACTION-DIRECTION':
			reverse = False
			if re.search('LEFT-TO-RIGHT',line[1]):
				reversible = False
			if re.search('RIGHT-TO-LEFT',line[1]):
				reversible = False
				reverse = True

		if line[0] == 'ENZYMATIC-REACTION':
			enames.append(line[1]) 

		if line[0] == 'SYNONYMS':
			commname.append(line[1])

		if line[0] == 'TYPES':
			typ.append(line[1])

		if line[0] == 'IN-PATHWAY':
			inpath.append(line[1])

		if line[0] == 'CANNOT-BALANCE?' and line[1] in ['T','NIL']:
			nobalance = True



	if reverse == False:
		for i in lines:
			line = i.strip().split(' - ')
			if line[0] == 'LEFT':
				innames.append(line[1].strip('|'))
				slst.append(1)
				tlst.append('CCO-IN')

			if line[0] == 'RIGHT':
				outnames.append(line[1].strip('|'))
				slst.append(-1)
				tlst.append('CCO-IN')

			if line[0] == '^COEFFICIENT':
				if 'n' in line[1]:
					line[1] = '0'
				if 'x' in line[1]:
					line[1] = '1'
				slst[len(slst)-1] = slst[len(slst)-1]*int(line[1])

			if line[0] == '^COMPARTMENT':
				tlst[len(tlst)-1] = line[1]
			


	if reverse == True:
		for i in lines:
			line = i.strip().split(' - ')
			if line[0] == 'RIGHT':
				innames.append(line[1].strip('|'))
				slst.append(1)
				tlst.append('CCO-IN')

			if line[0] == 'LEFT':
				outnames.append(line[1].strip('|'))
				slst.append(-1)
				tlst.append('CCO-IN')

			if line[0] == '^COEFFICIENT':
				if 'n' in line[1]:
					line[1] = '0'
				if 'x' in line[1]:
					line[1] = '1'
				slst[len(slst)-1] = slst[len(slst)-1]*int(line[1])

			if line[0] == '^COMPARTMENT':
				tlst[len(tlst)-1] = line[1]


	r = EC.Reaction([rname, commname, typ, innames, outnames, slst, tlst, nobalance, reversible, enames, inpath, dic])
	dic.appendr(r)


def estring(string,dic):
	""" splits and stores the strings coming from the regulation file """

	name = ''
	commname = ''
	regulator = ''
	renzname = ''
	mechanism = 'x'
	activation = None
	found = False
	cs = []
	rct = None
	lines = string.strip().split('\n')
	for i in lines:
		line = i.strip().split(' - ')
		if line[0] == 'TYPES':
			if line[1] == 'Regulation-of-Enzyme-Activity':
				found = True
				break

	if found == True:
		for i in lines:
			line = i.strip().split(' - ')
			if line[0] == 'REGULATED-ENTITY':
				renzname = line[1]
				rct = dic.renzsearch(renzname)
				found = False

		if found == True:
			return None

		for i in lines:
			line = i.strip().split(' - ')
			if line[0] == 'UNIQUE-ID':
				name = line[1]

		for i in lines:
			line = i.strip().split(' - ')
			if line[0] == 'COMMON-NAME':
				commname = line[1]


			if line[0] == 'MECHANISM' and mechanism == 'x':
				if line[1] == 'COMPETITIVE':
					mechanism = 'c'
				elif line[1] == 'ALLOSTERIC':
					mechanism = 'a'
				else:
					mechanism = 'n'
				
			if line[0] == 'MODE' and activation == None:
				if line[1] == '-':
					activation = False
				elif line[1] == '+':
					activation = True
				else:
					activation = None

			if line[0] == 'REGULATOR' and regulator == '':
				g = dic.gsearch(line[1])
				cs.append(g)
				if g == None:
					return None

		e = EC.Regulation(name, commname, cs, mechanism, activation, rct)
		rct.elst.append(e)
		dic.appende(e)


def pstring(temp, firstiteration, dic):
	""" splits and stores the strings coming from the pathways file """

	#rctlst = ''

	lines = temp.split('\n')
	for i in lines:
		l = i.strip().split(' - ')
		if l[0] == 'UNIQUE-ID':
			name = l[1]
			p = dic.psearch(name)
			if p == None:
				p = EC.Pathway(name)
				dic.appendp(p)
			if p.finished == True:
				return True
		if l[0] == 'COMMON-NAME' and p.commname == '':
			p.commname = l[1]
			
		if l[0] == 'REACTION-LAYOUT':	
			inel = []
			outel = []
			rct = re.search('\((.*) \(:LEFT-PRIMARIES(.*)\) \(:DIRECTION (.*)\) \(:RIGHT-PRIMARIES(.*)\)\)',l[1])
			if rct.group(2)=='' and rct.group(4)=='':
				subpname = rct.group(1)
				subp = dic.psearch(subpname)
				if subp == None:
					subp = EC.Pathway(subpname)
					dic.appendp(subp)
				if firstiteration:
					p.todo.append(subp)
					p.splst.append(subp)
				elif subp.finished == True and subp in p.todo:
					for i in range(len(subp.rlst)):
						p.appendr(subp.rlst[i],subp.primlst[i])
					p.todo.remove(subp)
					if len(p.todo) == 0:
					#the last subp has just been added: declare the path as finished and
					#add the pathway to all the plst of all compounds and to reactions involved
						p.finished = True
				else:
					continue				
			elif firstiteration and rct.group(3) == ':L2R':
				inel = dic.glst([i.strip('|') for i in rct.group(2).strip().split(' ')])
				outel = dic.glst([i.strip('|') for i in rct.group(4).strip().split(' ')])
				r = dic.rsearch(rct.group(1).strip('|'))
				if r == None:
					print('Fatal error during decoding of pathway '+p.name+': the reaction '+rct.group(1)+' has not been found in the database')
				p.appendr(r,(inel,outel))
				appinpath(inel,p)
				appinpath(outel,p)

			elif firstiteration and rct.group(3) == ':R2L':
				inel = dic.glst([i.strip('|') for i in rct.group(4).strip().split(' ')])
				outel = dic.glst([i.strip('|') for i in rct.group(2).strip().split(' ')])
				r = dic.rsearch(rct.group(1).strip('|'))
				if r == None:
					print('Fatal error during decoding of pathway '+p.name+': the reaction '+rct.group(1)+' has not been found in the database')
				p.appendr(r,(inel,outel))
				appinpath(inel,p)
				appinpath(outel,p)

			elif firstiteration:
				print 'broken pathway: '
				print temp
				p.todo = []
	
	if len(p.todo) == 0:
		p.finished = True
		return True
	else:
		return False
	
def mstring(line,colors,dic):
	""" 
	calculates all the possible color configurations on the edges contained in line,
	creates the corresponding object and stores them as motifs in the dic attribute mlst.
	"""
	lst = line.strip().split(' ')
	elst = [(int(lst[i]),int(lst[i+1])) for i in range(2,len(lst),2)]
	g = Graph(int(lst[0]),directed=True)
	g.vs['name'] = range(len(g.vs[:]))
	g.add_edges(elst)
	clst = colorlst(elst,colors)
	for c in clst:
		print c
		g.es['color'] = c
		m = EC.Motif(g.copy())
		dic.appendm(m)
		
def colorlst(elst,colors):
	""" gives a list of possible color configurations on the edges elst"""
	clst = [list(p) for p in itertools.product(colors,repeat=len(elst))]
	return clst

def appinpath(clst, p):
	for c in clst:
		c.inpath.append(p.name)
	
def ercheck(temp, dic):
	renzname = ''
	qname = ''
	rname = ''
	ename = ''
	lines = temp.split('\n')
	for i in lines:
		line = i.strip().split(' - ')
		if line[0] == 'UNIQUE-ID':
			renzname = line[1]
		if line[0] == 'ENZYME':
			qname = line[1]
			q = dic.qsearch(qname)
			if q == None:
				print('Registration error for '+renzname+':\t\t'+qname+' not found')
			elif renzname not in q.enz:
				print('Registration error for '+q.name+':\t\t'+str(q.enz)+' -- '+renzname)
		if line[0] == 'REACTION':
			rname = line[1]
			r = dic.rsearch(rname)
			if r == None:
				print('Registration error for '+renzname+':\t\t'+rname+' not found')
			elif renzname not in r.ename:
				print('Registration error for '+r.name+':\t\t'+str(r.ename)+' -- '+renzname)
		if line[0] == 'REGULATED-BY':
			ename = line[1]
			e = dic.esearch(ename)
			if e == None:
				print('Registration error for '+renzname+':\t\t'+ename+' not found')
			elif renzname not in e.rct.ename:
				print('Registration error for '+e.name+':\t\t'+str(e.rct.ename)+' -- '+renzname)


	


#interface
#add functions
#create graph

if __name__ == "__main__":
	organism = raw_input('Organism? (ecoli/yeast[custom]) \n')
	if organism == '':
		organism = 'yeast'
	dic = EC.ReactDic(organism)
	decode(organism,dic)
	
	with open('./results/logmeta.log','w') as out1:
		for r in dic.rlst:
			if len(r.elst) >= 2: 
				out1.write('\n'+r.name+'\n')	
				for e in r.elst:
					out1.write(str(e.activation)+'\t'+e.cs[0].name+'\n')
	ppcorrp = [0] * (len(dic.plst)*len(dic.plst))
	ppcorrn = [0] * (len(dic.plst)*len(dic.plst))
	for e in dic.elst:
		cs = e.cs[0]
		r = e.rct
		if e.activation == None:
			continue
		for (i,p1) in enumerate(dic.plst):
			for (j,p2) in enumerate(dic.plst):
				if cs in p1.sclst and r in p2.srlst and p1 not in p2.splst and p2 not in p1.splst and p1 != p2:
					if e.activation == True:
						ppcorrp[i*len(dic.plst)+j] += 1
					elif e.activation == False:
						ppcorrn[i*len(dic.plst)+j] += 1

	ppcorrt = ppcorrp
		
	with open('./results/logpath.log','w') as out2:
		for i in range(20):
			j = ppcorrp.index(max(ppcorrp))
			(m,n) = ((j-j%len(dic.plst))/len(dic.plst),j%len(dic.plst))
			(p1,p2) = (dic.plst[m],dic.plst[n])
			out2.write(str(ppcorrp[j])+'\t('+str(ppcorrp[n*len(dic.plst)+m])+')\t'+str([ppcorrn[j],ppcorrn[n*len(dic.plst)+m]])+'\n'+p1.name+'\t'+p2.name+'\n')
			ppcorrp[m*len(dic.plst)+n] = 0
			ppcorrp[n*len(dic.plst)+m] = 0


		out2.write('\n-------------------------------------------------------------------------------------\n')

		for i in range(20):
			j = ppcorrn.index(max(ppcorrn))
			(m,n) = ((j-j%len(dic.plst))/len(dic.plst),j%len(dic.plst))
			(p1,p2) = (dic.plst[m],dic.plst[n])
			out2.write(str(ppcorrn[j])+'\t('+str(ppcorrn[n*len(dic.plst)+m])+')\t'+str([ppcorrt[j],ppcorrt[n*len(dic.plst)+m]])+'\n'+p1.name+'\t'+p2.name+'\n')
			ppcorrn[m*len(dic.plst)+n] = 0
			ppcorrn[n*len(dic.plst)+m] = 0


	print('\nProgram executed correctly.')

