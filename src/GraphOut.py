#!/usr/bin/python

from igraph import *

def write_pathway_gml(gtot, elst, outfilename):
	""" 
	this function writes a gml file containing the grapical 
	representation of the graph gtot by highliting and spring-embedding 
	the nodes of the principal pathway and adding the rest of the nodes 
	of the pathway at a second time. 
	"""

	gprim = gtot.copy()
	gprim.delete_vertices([i for i,inprim in enumerate(gtot.vs['inprim'][:]) if inprim == 0])
	tlayout = gprim.layout_kamada_kawai()
	primlayout = gprim.layout_graphopt(niter=500,node_charge=0.07,node_mass=30,spring_constant=15,spring_length=35,max_sa_movement=5,seed=tlayout)

	with open(outfilename,'w') as outf:
		outf.write('Version 1\ngraph\n[\n\tdirected 1\n')
		# write nodes and their position
		for i,n in enumerate(gtot.vs['name'][:]):
			primid = getid(gprim.vs['name'][:],n)
			if primid != None:
				outf.write('\tnode\n\t[\n\t\tid '+str(i)+'\n\t\tname "'+n+'"\n\t\tinpath 1\n\t\tpnamelst "'+gtot.vs['pnamelst'][i]+'"\n\t\tgraphics\n\t\t[\n\t\t\tx '+str(primlayout[primid][0])+'\n\t\t\ty '+str(primlayout[primid][1])+'\n')
				if 'RXN' in n:
					outf.write('\t\t\ttype "rectangle"\n\t\t\tw 40\n\t\t\th 40\n\t\t\tfill "#FF0000"\n\t\t]\n\t]\n')
				else:
					outf.write('\t\t\ttype "circle"\n\t\t\tw 40\n\t\t\th 40\n\t\t\tfill "#FF0000"\n\t\t]\n\t]\n')
			else:
				outf.write('\tnode\n\t[\n\t\tid '+str(i)+'\n\t\tname "'+n+'"\n\t\tinpath 1\n\t\tpnamelst "'+gtot.vs['pnamelst'][i]+'"\n\t\tgraphics\n\t\t[\n\t\t\tx '+str(0)+'\n\t\t\ty '+str(0)+'\n')
				if 'RXN' in n:
					outf.write('\t\t\ttype "rectangle"\n\t\t\tw 40\n\t\t\th 40\n\t\t\tfill "#0000FF"\n\t\t]\n\t]\n')
				else:
					outf.write('\t\t\ttype "circle"\n\t\t\tw 40\n\t\t\th 40\n\t\t\tfill "#0000FF"\n\t\t]\n\t]\n')
		for e in gtot.get_edgelist():
			outf.write('\tedge\n\t[\n\t\tsource '+str(e[0])+'\n\t\ttarget '+str(e[1])+'\n\t\tfill "#000000"\n\t\ttype "reaction"\n\t\ttargetArrow "delta"\n\t]\n')
		for e in elst:
			cids = [getid(gtot.vs['name'][:],c.name) for c in e.cs if getid(gtot.vs['name'][:],c.name) != None]
 			#Forward reaction
			rid = getid(gtot.vs['name'][:],e.rct.name)
			if rid != None:
				for cid in cids:
					if e.activation == True:
						outf.write('\tedge\n\t[\n\t\tsource '+str(cid)+'\n\t\ttarget '+str(rid)+'\n\t\ttype "positive-regulation"\n\t\tgraphics\n\t\t[\n\t\t\tsmoothBends 1\n\t\t\tfill "#00FFFF"\n\t\t\ttargetArrow "delta"\n\t\t]\n\t]\n')
					if e.activation == False:
						outf.write('\tedge\n\t[\n\t\tsource '+str(cid)+'\n\t\ttarget '+str(rid)+'\n\t\ttype "negative-regulation"\n\t\tgraphics\n\t\t[\n\t\t\tsmoothBends 1\n\t\t\tfill "#FF0000"\n\t\t\ttargetArrow "delta"\n\t\t]\n\t]\n')

 			#Reverse reaction
			rid = getid(gtot.vs['name'][:],e.rct.name+'_r')
			if rid != None:
				for cid in cids:
					if e.activation == True:
						outf.write('\tedge\n\t[\n\t\tsource '+str(cid)+'\n\t\ttarget '+str(rid)+'\n\t\ttype "positive-regulation"\n\t\tgraphics\n\t\t[\n\t\t\tsmoothBends 1\n\t\t\tfill "#00FFFF"\n\t\t\ttargetArrow "delta"\n\t\t]\n\t]\n')
					if e.activation == False:
						outf.write('\tedge\n\t[\n\t\tsource '+str(cid)+'\n\t\ttarget '+str(rid)+'\n\t\ttype "negative-regulation"\n\t\tgraphics\n\t\t[\n\t\t\tsmoothBends 1\n\t\t\tfill "#FF0000"\n\t\t\ttargetArrow "delta"\n\t\t]\n\t]\n')



		outf.write(']')

def write_hub_gml(gtot, elst, outfilename):
	""" 
	this function writes a gml file containing the grapical 
	representation of the graph ghub by highliting and spring-embedding 
	the hub nodes and adding the activation/repression intractions at a 
	later time. 
	"""

	tlayout = gtot.layout_kamada_kawai()
	layout = gtot.layout_graphopt(niter=500,node_charge=0.07,node_mass=30,spring_constant=15,spring_length=35,max_sa_movement=5,seed=tlayout)

	with open(outfilename,'w') as outf:
		outf.write('Version 1\ngraph\n[\n\tdirected 1\n')
		# write nodes and their position
		for i,n in enumerate(gtot.vs['name'][:]):
			outf.write('\tnode\n\t[\n\t\tid '+str(i)+'\n\t\tname "'+n+'"\n\t\tgraphics\n\t\t[\n\t\t\tx '+str(layout[i][0])+'\n\t\t\ty '+str(layout[i][1])+'\n')
			if 'RXN' in n or 'PWY' in n:
				outf.write('\t\t\ttype "rectangle"\n\t\t\tw 40\n\t\t\th 40\n\t\t\tfill "#FF0000"\n\t\t]\n\t]\n')
			else:
				outf.write('\t\t\ttype "circle"\n\t\t\tw 40\n\t\t\th 40\n\t\t\tfill "#0000FF"\n\t\t]\n\t]\n')

		for e in gtot.get_edgelist():
			outf.write('\tedge\n\t[\n\t\tsource '+str(e[0])+'\n\t\ttarget '+str(e[1])+'\n\t\tfill "#000000"\n\t\ttype "reaction"\n\t\ttargetArrow "delta"\n\t]\n')
		for e in elst:
			cids = [getid(gtot.vs['name'][:],c.name) for c in e.cs if getid(gtot.vs['name'][:],c.name) != None]
 			#Forward reaction
			rid = getid(gtot.vs['name'][:],e.rct.name)
			if rid != None:
				for cid in cids:
					if e.activation == True:
						outf.write('\tedge\n\t[\n\t\tsource '+str(cid)+'\n\t\ttarget '+str(rid)+'\n\t\ttype "positive-regulation"\n\t\tgraphics\n\t\t[\n\t\t\tsmoothBends 1\n\t\t\tfill "#00FFFF"\n\t\t\ttargetArrow "delta"\n\t\t]\n\t]\n')
					if e.activation == False:
						outf.write('\tedge\n\t[\n\t\tsource '+str(cid)+'\n\t\ttarget '+str(rid)+'\n\t\ttype "negative-regulation"\n\t\tgraphics\n\t\t[\n\t\t\tsmoothBends 1\n\t\t\tfill "#FF0000"\n\t\t\ttargetArrow "delta"\n\t\t]\n\t]\n')

		outf.write(']')


			



def getid(lst,string):
	for i,x in enumerate(lst):
		if x == string:
			return i
	return None

if __name__ == "__main__":
	g1 = Graph(0)
	g1.add_vertices(['a','b','c'])
	g1.add_edges([('a','b'),('b','c'),('c','a')])
	g1.vs['pnamelst'] = ['x','q','x|q|r']
	g1.vs['inprim'] = [1,0,1]
	write_pathway_gml(g1,[],'./tests/GraphOutTest.gml')
