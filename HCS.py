
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
	#Sink for Edward-Karp Algorithm
	end = edge_list[len(edge_list)-1][0]
	if start !=end:
		edge_set = minimum_st_edge_cut(g,start,end)
	else:
		edge_set = []
	a = len(edge_set)
	b = len(edge_list)
	
	return 1
	



""" HCS Algorithm for Recursive Calls """
def HCS(g):
	return 1


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

 



""" Function to extract each connected component of the Graph (in form of edges) & pass it onto Minimum Cut"""
def main_clustering(component, main_graph):
	#Forming a temporary Graph Structure for the connected component
	temp_graph = nx.Graph()
	temp_graph = main_graph.copy()
	for item in temp_graph.edges():
		#First Node
		source = item[0]
		#Second Node
		target = item[1]
		if source not in component and target not in component:
			#Removal this edge from the Data Structure
			temp_graph.remove_edge(source,target)


	#End of removal of redundant edges
	#print len(temp_graph.edges())

	#Find the edges involved in minimum cut
	cut_edges = find_minimum_cut(temp_graph)
	





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
























