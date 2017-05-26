
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Application of Highly Connected Sub-graph  [HCS] algorithm on the Preliminary Dataset 
   -> Treatment of the Data Set as Undirected Graph -- Just Interactions taken into account """


""" Highly Connected Subgraph Algorithm : 
	  1. Find minimum cut 
	  2. Obtain Sub-graphs
	  3. Check which one of the subgraphs are not highly Connected
	  4. Return to Step 1 for those subgraphs                                                 """


#Libraries 
import json
import pandas as pd
from py2neo import Node, Relationship, Graph
import re
import string
import networkx as nx
from networkx.algorithms.connectivity import minimum_st_edge_cut
import matplotlib.pyplot as plt

""" Max Flow algorithm for Minimum Cut """
#Function to find the Minimum Cut of a Graph
def find_minimum_cut(g):

	#List of edges in the component
	edge_list = g.edges()
	#Source for Edward-Karp Algorithm
	start = edge_list[0][0]
	#Sink Initialization
	end = edge_list[0][0]
	#Sink for Edward-Karp Algorithm -- Sink Reassignment s.t : Source != Sink
	for item in edge_list:
		if item[0] != start:
			end = item[0]
			break
		elif item[1]!=start:
			end = item[1]
			break
		else:
			continue
			
			

	if start != end:
		#List of edges comprising of the Cut Set
		edge_set = minimum_st_edge_cut(g,start,end)
	else:
		edge_set = []


	#List of edges in cut-set
	cut_set = list(edge_set)

	return cut_set

	
	
	


#Function to check the presence of any self loops -- If present, they have to be removed
def self_loop_presence(g):
	#Looping through edge
	for item in g.edges():
		if item[0] == item[1]:
			#At this stage remove the self loop from the Graph Data Structure
			g.remove_edge(item[0],item[1])
			#print "Found"
			

	#Return Graph Data Structure without self-loops
	return g


clusters = []

#Function for Highly Connected Sub-graph Formation
def HCS(connected_graph):
	#Initialization for connected sub-graphs
	temp_connected_graph = nx.Graph()
	temp_connected_graph = connected_graph

	#Find the Minimum Cut-edges
	cut_edges = find_minimum_cut(temp_connected_graph)

	total =  len(temp_connected_graph.edges())

	number_of_nodes = len(temp_connected_graph.nodes())

	edge_connectivity = len(cut_edges)

	""" For a Graph to be HCS -- edge connectivity > number_of_nodes/2  """

	#Condition True : HCS is valid
	if edge_connectivity > (number_of_nodes/2): 
		#Store the nodes in a temporary list
		temp_list = []
		temp_list.append(temp_connected_graph.nodes)
		clusters.append(temp_list)

    
    #Divide the two components and pass into HCS
	else:
		#Removal of the edges from the Graph to form separate connected components
		for item in cut_edges:
			temp_connected_graph.remove_edge(item[0],item[1])


		connected_components = list(nx.connected_components(temp_connected_graph))

		#At this stage we will get Two Connected Components -- Compute & Store them separately
		first_graph = nx.Graph()
		second_graph = nx.Graph()

		first_graph = temp_connected_graph.copy()
		second_graph = temp_connected_graph.copy()

		""" Separation of components """
		#First Component
		first_component = connected_components[0]
		#Second Component
		second_component = connected_components[1]

		#Removal of Nodes - For First Component
		for item in first_graph.nodes():
			if item not in first_component:
				first_graph.remove_node(item)

		#Removal of Edges - For First Component
		for item in first_graph.edges():
			#First Node
			source = item[0]
			#Second Node
			target = item[1]
			#Nodes corresponding to the current edge does not belong to the first component
			if source not in first_component and target not in first_component:
				first_graph.remove_edge(source,target)


		#Removal of Nodes - For Second Component
		for item in second_graph.nodes():
			if item not in second_component:
				second_graph.remove_node(item)


		#Removal of Edges - For Second Component
		for item in second_graph.edges():
			#First Node
			source = item[0]
			#Second Node
			target = item[1]
			#Nodes corresponding to the current edge does not belong to the second component
			if source not in second_component and target not in second_component:
				second_graph.remove_edge(source,target)



		""" Testing Criterion : total = len(cut_edges) + len(first_graph.edges()) + len(second_graph.edges()) """

		#if total != (len(cut_edges) + len(first_graph.edges()) + len(second_graph.edges())):
		#	print "Error"

        #print len(list(nx.connected_components(second_graph)))
		#print len(list(nx.connected_components(second_graph)))


		HCS(first_graph)
		HCS(second_graph)



		







 



""" Function to extract each connected component of the Graph (in form of edges) & pass it onto Minimum Cut"""
def main_clustering(component, main_graph):
	#Forming a temporary Graph Structure for the connected component
	temp_graph = nx.Graph()
	temp_graph = main_graph.copy()


	""" According to Network Graph Data-structure -- it is essential to remove both nodes and edges to modify a graph """

	#Removal of the nodes from temp_graph which are present not in component
	for item in temp_graph.nodes():
		if item not in component:
			temp_graph.remove_node(item)

	
	#Removal of edges from the temp_graph which are not present in component
	for item in temp_graph.edges():
		#First Node
		source = item[0]
		#Second Node
		target = item[1]
		if source not in component and target not in component:
			#Removal this edge from the Data Structure
			temp_graph.remove_edge(source,target)



	#End of removal of redundant edges
	

	#Find the edges involved in minimum cut -- "They have to be removed to obtain two components"
	cut_edges = find_minimum_cut(temp_graph)

	

	""" Tests : Number of Connected Components == 1 """
	#print len(list(nx.connected_components(temp_graph)))	

	""" At this stage we have a connected component & a list of cut-edges on whose removal the graph will become disconnected """

	#print cut_edges

	""" Tests : Number of Connected Components == 2"""
	#for item in cut_edges:
	#	temp_graph.remove_edge(item[0],item[1])

	#print len(list(nx.connected_components(temp_graph)))

	
	#Calling HCS function : @parameters => #temp_graph : Connected Graph
	HCS(temp_graph)









	





#Loading gene_interactions JSON file into a variable 
with open('JSON rows/gene_interactions.json') as json_data:
	interactions = json.load(json_data)


#Information about Graph Connectivity
graph_info = interactions["results"]

source = []
target = []

#Creating NetworkX instance
graph = nx.Graph()

#Extracting the edge relationships
for edge in graph_info:
	source.append(edge[0])
	target.append(edge[2])

	#Adding the edge in NetworkX
	graph.add_edge(edge[0],edge[2])



#print len(graph.nodes())
#print len(graph.edges())
#Obtaining the Graph Data Structure after removal of Self-Loops
graph = self_loop_presence(graph)


""" Commented Code - For drawing the Graphs """
#Drawing the Graph in MatPlotLib   --- Takes a lot of time to load
#nx.draw(graph)
#nx.draw_random(graph)
#nx.draw_circular(graph)
#nx.draw_spectral(graph)
#plt.show()



#List of Connected Components
connected_components =  list(nx.connected_components(graph))

""" All the individual connected components have to be analysed -- For the Preliminary Data set --96 components with two Major Components """


#Initial Number of Clusters :
initial_cluster = len(connected_components)

#print initial_cluster

#print connected_components[2]

connected_components_as_subgraph = nx.connected_component_subgraphs(graph,copy=True)
g = list(connected_components_as_subgraph)

#Looping to pass each component into main_clustering to obtain sub_clusters
for component in g:
	#Calls main_clustering each time with connected components : list(component) - (List of nodes in that component) & Original Graph Structure 
	main_clustering(list(component),graph)
	#print len(list(component))
	break
	








