#!/usr/bin/python

from igraph import Graph

class ReactDic(object):
	
	# graphname = name of the graph
	# clst = list of compund objects ordered by their ID	
	# rlsts = list of reaction objects listed by their IDs
	# elst = list of enzyme objects listed by their IDs
	# plst = list of pathways listed by their IDs
	# mlst = list of motifs

	def __init__(self,name):
		self.graphname = name
		self.clst = []
		self.rlst = []
		self.elst = []
		self.plst = []
		self.klst = []
		self.qlst = []
		self.mlst = []
		self.radjlst = []
		self.dgraph = Graph(0,directed = True)

	# gives back the name of the graph
	def getname(self):
		return self.graphname



	# the following three functions count the number of elements in the respective lst
	def ccount(self):
		return len(self.clst)

	def rcount(self):
		return len(self.rlst)

	def ecount(self):
		return len(self.elst)

	def pcount(self):
		return len(self.plst)



	# the following three functions add an element to the respective lst and give back its ID
	def appendc(self,c):
		self.clst.append(c)
		return len(self.clst)-1

	def appendr(self,r):
		self.rlst.append(r)
		return len(self.rlst)-1

	def appende(self,e):
		self.elst.append(e)
		return len(self.elst)-1

	def appendp(self,p):
		self.plst.append(p)
		return len(self.plst)-1

	def appendk(self,k):
		self.klst.append(k)
		return len(self.klst)-1

	def appendq(self,q):
		self.qlst.append(q)
		return len(self.qlst)-1
	
	def appendm(self,m):
		self.mlst.append(m)
		return len(self.mlst)-1



	# search in the dictionary dic for name and give back the corresponding instance. Returns None if no instance has been found.
	def rsearch(self,rname):
		for r in self.rlst:	
			if rname.strip('|') == r.name:
				return r

	def csearch(self,cname):
		for c in self.clst:	
			if cname == c.name:
				return c

	def esearch(self,ename):
		for e in self.elst:	
			if ename == e.name:
				return e

	def psearch(self,pname):
		for p in self.plst:
			#print [p.__dict__ for p in self.plst]
			if pname == p.name:
				return p

	def ksearch(self,kname):
		for k in self.klst:	
			if kname == k.name:
				return k

	def qsearch(self,qname):
		for q in self.qlst:	
			if qname == q.name:
				return q

	def renzsearch(self,ename):
		for r in self.rlst:
			if ename in r.ename:
				return r
			
	def msearch(self,name):
		for m in self.mlst:
			if name in m.name:
				return m

	def gsearch(self,gname):
		gn = gname.strip('|')
		g = self.csearch(gn)
		if g == None:
			g = self.qsearch(gn)
			if g == None:
				g = self.ksearch(gn)
				if g == None:
					print('ERROR IN GSEARCH:    '+gname+' has not been found in the database.')
					g = Compound(gn, '', 'NOT FOUND', '', '', [])
					self.appendc(g)

					
		return g

	def glst(self,gnames):
		glst = []
		for gn in gnames:
			glst.append(self.gsearch(gn))
		return glst



	def createadjlst(self):
		""" initializes self.radjlst: the reaction adjacency list """
		self.radjlst = []
		for r1 in self.rlst:
			for r2 in self.rlst:
				for o in r1.outnames:
					for i in r2.innames:
						if o == i:
							self.radjlst.append([r1,r2])

	def children(self,r):
		""" returns the OUT r-neighbours of a reaction r """
		children = []
		for radj in self.radjlst:
			if radj[0] == r:
				children.append(radj[1])
		return children

	def kfillup(self):
		for k in self.klst:
			lst = []
			stack = k.children
			while stack != []:
				lst.extend(stack)
				stack2 = []
				for k2 in stack:
					stack2.extend(k2.children)
				stack = stack2[:]
			k.subclasses = lst

	def qfillup(self):
		"""substitutes the relations of proteins with the actual objects in self.qlst"""
		for q in self.qlst:
			compt = []
			compoft = []
			for qn in q.comp:
				compt.append(self.qsearch(qn))
			for qn in q.compof:
				compoft.append(self.qsearch(qn))
			q.comp = compt[:]
			q.compof = compoft[:]
	

	def kcs(self):
		for c in self.clst:
			for k in c.types:
				k.cs.append(c)

	def pflat(self):
		"""once the decoding of pathways has been carried out this function fills the single reactions and single primary elements of all the pathways"""
		for p in self.plst:
			p.srlst = list(set(p.rlst))
			for r in p.srlst:
				#print r
				#print r.__dict__
				r.plst.append(p)
			for prim in p.primlst:
				for i in range(2):
					for c in prim[i]:
						if c not in p.sclst:
							p.sclst.append(c)
							c.cplst.append(p)
				
			

	def create_dgraph(self):
		""" initializes self.graph: the total directed graph of reactions. It ereases any preceding graph. """
		self.dgraph = Graph(0,directed = True)
		cnlst = []
		rnlst = []
		rlinks = []
		for rct in self.rlst:
			rnlst.append(rct.name)	
			for i in rct.innames:
				cnlst.append(i)
				rlinks.append((i,rct.name))
			for o in rct.outnames:
				cnlst.append(o)
				rlinks.append((rct.name,o))
			if rct.reversible == True:
				rnlst.append(rct.name+'_r')
				for o in rct.innames:
					rlinks.append((rct.name+'_r',o))
				for i in rct.outnames:
					rlinks.append((i,rct.name+'_r'))
		cnlst = list(set(cnlst))
		cinpathlst = [self.gsearch(cn).inpath if self.gsearch(cn) != None else [] for cn in cnlst]
		rnlst = list(set(rnlst))
		rinpathlst = [self.rsearch(rn).inpath if '_r' not in rn else self.rsearch(rn[0:len(rn)-2]) for rn in rnlst]
		cnlst.extend(rnlst)
		cinpathlst.extend(rinpathlst)
		self.dgraph.add_vertices(cnlst)
		self.dgraph.add_edges(rlinks)
		self.dgraph.vs['inpath'] = cinpathlst
		#print rlinks[len(rlinks)-3000:len(rlinks)-1]
		print self.dgraph.summary()

	
	def create_pgraphs(self,filterlst):
		""" create the graphs of all the pathways and initialize the corresponding variable """
		for p in self.plst:
			p.dependence(self)
		print '\nprimgraph'
		for p in self.plst:
			p.create_primgraph(filterlst)
		print '\ntotgraph\n'
		for p in self.plst:
			p.create_totgraph(filterlst)

	def create_filterlst(self,threshold):
		""" create a filtering list nlst containing all the vertex in dgraph that have a degree higher than threshold """

		lst = []
		nlst = []
		for v in range(self.dgraph.vcount()):
			if self.dgraph.degree(v) > threshold:
				lst.append(v)
				nlst.append(self.dgraph.vs['name'][v])
		print 'NUMBER OF ELEMENTS TO FILTER: '+str(len(lst))
		print nlst
		return nlst

	def stoichmat(self):
		""" outputs the full stoichiometric matrix of the system """

		smat = []
		clst = []
		rlst = []
		rev = []

		for r in self.rlst:
			if r.nobalance == False:
				if r.name not in rlst:
					rlst.append(r.name)
					rev.append(r.reversible)
				else:
					print 'error in STOICHMAT: reaction\t',r.name,'\talready present!'
			for cname in union(r.innames[:],r.outnames[:]):
				if cname not in clst:
					clst.append(cname)

		for r in self.rlst:
			if r.nobalance == False:
				#create the stoichiometric vector
				rclst = r.innames[:]
				rclst.extend(r.outnames[:])
				sv = [r.slst[[i for i,b in enumerate(rclst) if b == cname][0]] if cname in rclst else 0 for cname in clst]
				smat.append(sv)

		return [smat,clst,rlst,rev]	
	
				


class Class(object):
	""" type object, defines a hierarchical structure in the compounds """

	def __init__(self,name):
		""" define the name of the class, its children (types that also belong to that class) and its common names """
		self.name = name
		self.children = []
		self.fathers = []
		self.commname = []
		self.subclasses = []
		self.gs = []
		self.cplst = []
		self.inpath = []

	def addchild(self,child):
		""" adds a child to the corresponding variable, one needs to search for the child in the dictionary first """
		self.children.append(child)

	def addfather(self,father):
		""" adds a child to the corresponding variable, one needs to search for the child in the dictionary first """
		self.fathers.append(father)

	def addg(self,g):	
		""" adds a general child to the corresponding variable, the child can be a compound or a protein (class object)"""
		self.gs.append(g)

	def setcommname(self, cmname):
		""" defines the variable self.commname """
		self.commname = cmname


class Compound(object):
	""" Compound class """
	
	def __init__(self, name, commname, typ, smiles, weight, reglst):
		""" 
		Initializes the variables of the object.
		Name is needed, all other variables will be initialized later if needed.
		"""
		self.name = name
		self.commname = commname
		self.typ = typ
		self.weitht = weight
		self.smiles = smiles
		self.reglst = reglst
		self.cplst = []
		self.inpath = []
		self.graph = Graph(0,directed = True)

	def create_cplst(self,dic):
		""" initializes the self.cplst variable: the list of pathways passing through c """
		for p in dic.plst:
			if self.name in p.primlstflat:
				self.cplst.append(p)

	def create_graph(self,filterlst):
		""" initializes the self.graph variable: the graph of pathways passing through c """
		for p in self.cplst:
			grapht = merge(self.graph,p.totgraph)
			self.graph = grapht

class Protein(object):
	""" protein class"""

	def __init__(self, name, commname, typ, comp, compof, reg, enz):
		"""
		Initializes the values of the object.
		Requires all the variables, that can be empty lists:
		name: UNIQUE ID
		commname: COMMON NAME
		typ: TYPE
		comps: COMPONENTS
		reg: REGULATES
		enz: CATALYZES
		"""
		self.name = name
		self.commname = commname
		self.typ = typ
		self.comp = comp
		self.compof = compof
		self.reg = reg
		self.enz = enz
		self.cplst = []
		self.inpath = []

class Reaction(object):
	""" reaction class""" 
	
	def __init__(self, lst):
		""" 
		Initializes the variables of the object.
		Name is needed, all other variables will be initialized later.
		"""
		[name, commname, typ, innames, outnames, slst, tlst, nobalance, reversible, enames, inpath, dic] = lst
		self.name = name
		self.commname = commname
		self.typ = typ
		self.innames = innames
		self.outnames = outnames
		self.reversible = reversible
		self.trans = False
		self.ename = enames
		self.inpath = inpath
		self.nobalance = nobalance
		self.slst = slst
		self.tlst = tlst
		self.istransport = (tlst == ['CCO-IN']*len(tlst))
		self.plst = []
		self.elst = []
		self.fblst = []
		self.inels = dic.glst(innames)
		self.outels = dic.glst(outnames)

class Regulation(object):
	""" Enzyme class """
	def __init__(self, name, commname, cs, competitive, activation, rct):
		""" 
		initialization of the enzyme object. 
		Needed is the Enzyme name. Addition of the catalyzed reaction and the activation compound will come later. 
		"""	
		self.name = name
		self.commname = commname
		self.cs = cs
		self.competitive = competitive  # 'c' for competitive, 'n' for noncompetitive and 'a' for allosteric 
		self.activation = activation # can have the values True, False or None
		self.rct = rct # associated reaction
		self.fbplst = [] # ???
		self.ffplst = [] # ???

	def fbplsts(self,plst1,plst2):
		self.fbplst = [plst1,plst2]

	def ffplsts(self,plst1,plst2):
		self.ffplst = [plst1,plst2]

	def stab(self):
		cnames = [c.name for c in self.cs]
		if intersection(cnames,self.rct.outnames) == [] and intersection(cnames,self.rct.innames) == []:
			return None
		elif  intersection(cnames,self.rct.innames) != [] and self.activation == False or intersection(cnames,self.rct.outnames) != [] and self.activation == True:
			return False
		elif  intersection(cnames,self.rct.innames) != [] and self.activation == True or intersection(cnames,self.rct.outnames) != [] and self.activation == False:
			return True
		else:
			return None

class Pathway(object):
	""" Pathway object """
	def __init__(self, name):
		""" 
		Initializes the variables of the object.
		Name is needed, all other variables will be initialized later if needed.
		"""
		self.name = name
		self.commname = ''
		self.rlst = [] #list of principal reactions, aligned with primlst
		self.elst = [] #list of enzymes catalyzing the pathway
		self.srlst = [] 
		self.rnlst = []
		self.primlst = [] #list of (inel,outel), where inel/outel are list of strings themselves
		self.sclst = []
		self.splst = []
		self.todo = []
		self.substrates = [] #counted with multiplicities
		self.byproducts = [] #counted with multiplicities
		self.finished = False
		self.primgraph = Graph(0,directed = True)
		self.totgraph = Graph(0,directed = True)
		
	def dependence(self, dic):
		""" initializes the variables substrates, byproducts and elst: the dependence of the pathway """
		for (j,prim) in enumerate(self.primlst):
			r = self.rlst[j]
			if prim[0][0] in r.inels:
				self.substrates.extend(r.inels)
				self.byproducts.extend(r.outels)
			else:
				self.substrates.extend(r.outels)
				self.byproducts.extend(r.inels)
			self.elst.extend(r.elst)
				
	def appendr(self,r,prim):
		""" append reaction and its primaries to the corresponding lists """
		self.rlst.append(r)
		self.primlst.append(prim)

	def create_primgraph(self, filterlst):
		""" initializes the self.primgraph variable: the primary component graph of the pathway p """
		for (j,prim) in enumerate(self.primlst):
			r = self.rlst[j]
			vsearch(self.primgraph,r.name)
			for i in prim[0]:
				vsearch(self.primgraph,i.name)
				self.primgraph.add_edges([(i.name,r.name)])
			for o in prim[1]:
				vsearch(self.primgraph,o.name)
				self.primgraph.add_edges([(r.name,o.name)])

		if self.primgraph.vcount() != 0:
			vn = [v for v in self.primgraph.vs['name'][:] if v in filterlst]
			if vn != []:
				self.primgraph.delete_vertices(vn)


	def create_totgraph(self, filterlst):
		""" initializes the self.totgraph variable: the total graph of the pathway p """
		for r in self.srlst:
			vsearch(self.totgraph,r.name)
			for i in r.innames:
				vsearch(self.totgraph,i)
				self.totgraph.add_edges([(i,r.name)])
			for o in r.outnames:
				vsearch(self.totgraph,o)
				self.totgraph.add_edges([(r.name,o)])
			if r.reversible == True:
				vsearch(self.totgraph,r.name+'_r')
				for i in r.innames:
					self.totgraph.add_edges([(r.name+'_r',i)])
				for o in r.outnames:
					self.totgraph.add_edges([(o,r.name+'_r')])


		if self.totgraph.vcount() != 0:
			for vn in filterlst:
				if vn in self.totgraph.vs['name'][:]:
					self.totgraph.delete_vertices([vn])

		inprimlst = [1 if name in self.primgraph.vs['name'][:] or name[0:len(name)-2] in self.primgraph.vs['name'][:] else 0 for name in self.totgraph.vs['name']]
		self.totgraph.vs['inprim'] = inprimlst
		pnamelst = [self.name if int(i) == 1 else '' for i in inprimlst]
		self.totgraph.vs['pnamelst'] = pnamelst
		
		
class Network(object):
	"""
	The graphs associated with the network of chemical reactions
	This class is constructed in order to speed up the computational time needed to work with graphs
	In order to use it correctly, the attributes of every object of this class NEED to be updated after every modification!!!
	Examples
	- Directed graph
	- Flux Coupling Graph
	- ...
	"""	
	def __init__(self,name):
		"""
		Attributes of the class:
		- name = name
		- graph = element of the graph class from the package igraph. Inherits all the methods of the class
		- vlst = vertex ids list
		- vnlst = vertex names list (if present)
		- elst = edge list in the form (id1,id2)
		- eclst = edge color list (if present)
		"""
		self.name = name
		self.graph = Graph(0,directed=True)
		self.vlst = []
		self.vnlst = []
		self.elst = []
		self.eclst = []
		
	def update(self,graph):
		"""
		update the information about the element "self" of the Network class after elaboration 
		"""
		self.graph = graph
		self.vlst = range(len(graph.vs[:]))
		if 'name' in graph.vertex_attributes():
			self.vnlst = graph.vs['name'][:]
		self.elst = graph.get_edgelist()
		if 'color' in graph.edge_attributes():
			self.eclst = graph.es['color'][:]
				
class Motif(object):
	"""
	a class containing the graphs of the network motifs
	"""
	def __init__(self,graph):
		"""
		Attributes of the class:
		- graph = element of the graph class from the package igraph. Inherits all the methods of the class
		- vlst = vertex ids list
		- elst = edge list in the form (id1,id2)
		- eclst = edge color list (if present)
		"""
		self.graph = graph
		self.vlst = range(len(graph.vs[:]))
		self.elst = graph.get_edgelist()
		if 'color' in graph.edge_attributes():
			self.eclst = graph.es['color'][:]
		else:
			self.graph.es['color'] = [0]*len(self.elst)
		self.name = ','.join(['-'.join([str(i[0]),str(i[1])]) for i in self.elst])
		self.olst = []
		
	def store_olst(self,gname,lst):
		"""
		store the occurrence list of motif self in the network gname
		"""
		self.olst.append((gname,lst))

	def graphname(self,graph):
		return graph
		
def vsearch(g,vn):
	""" search for a vertex with name vn in the graph g. If the vertex is not in the graph then add it to the graph and return its ID. """
	try:
		ver = [i for i,v in enumerate(g.vs['name']) if v == vn]
	except KeyError:
		ver = []
	if ver == []:
		g.add_vertices([vn])
		return g.vcount()
	else:
		return ver[0]

def merge(g1,g2):
	""" merges graph g1 and graph g2 into the output graph"""
	g3nslst = list(set(g1.vs['name'][:]) | set(g2.vs['name'][:])) 
	g3 = Graph(0,directed=True)
	g3.add_vertices(g3nslst)
	g3elst = []
	for e in g1.get_edgelist():
		g3elst.append((g1.vs['name'][e[0]],g1.vs['name'][e[1]]))
	for e in g2.get_edgelist():
		g3elst.append((g2.vs['name'][e[0]],g2.vs['name'][e[1]]))
	g3.add_edges(g3elst)
	g3.simplify()
	#add attributes
	g1primlst = [vn for i,vn in enumerate(g1.vs['name'][:]) if int(g1.vs['inprim'][i]) == 1]
	g2primlst = [vn for i,vn in enumerate(g2.vs['name'][:]) if int(g2.vs['inprim'][i]) == 1]
	g3prim = [1 if vn in g1primlst or vn in g2primlst else 0 for vn in g3.vs['name'][:]]
	g3pnamelst = [[] for i in range(len(g3.vs['name'][:]))]
	for i,vn1 in enumerate(g3.vs['name'][:]):
		for j,vn2 in enumerate(g1.vs['name'][:]):
			if vn1 == vn2:
				g3pnamelst[i].extend(g1.vs['pnamelst'][j].strip().split('|'))
		for j,vn2 in enumerate(g2.vs['name'][:]):
			if vn1 == vn2:
				g3pnamelst[i].extend(g2.vs['pnamelst'][j].strip().split('|'))
	g3.vs['pnamelst'] = ['|'.join(map(str,list(set(inp)))) if inp != [] else '' for inp in g3pnamelst]
	#print g1.vs['pnamelst'][:]	
	#print g3.vs['name'][:]
	g3.vs['inprim'] = g3prim
	return g3

def intersection(a,b):
	lst = []
	for i in a:
		if i in b:
			lst.append(i)
	return lst

def union(a,b):
	lst = a[:]
	lst.extend(b[:])
	return lst


		


"""
 This piece was in create_pringraph in the pathway class
				for prim in self.primlst:
					if prim[0][0] in r.innames and prim[1][0] in r.outnames:
						found = True
						vsearch(self.graph,r.reactname)
						for i in r.innames:
							vsearch(self.graph,i)
							self.graph.add_edges([(i,r.reactname)])
						for o in r.outnames:
							vsearch(self.graph,o)
							self.graph.add_edges([(r.reactname,o)])
					elif prim[0][0] in r.outnames and prim[1][0] in r.innames:
						found = True
						vsearch(self.graph,r.reactname+'_r')
						for i in r.innames:
							vsearch(self.graph,i)
							self.graph.add_edges([(r.reactname+'_r',i)])
						for o in r.outnames:
							vsearch(self.graph,o)
							self.graph.add_edges([(o,r.reactname+'_r')])
"""

