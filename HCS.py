
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
import matplotlib.pyplot as plt

""" Karger's algorithm for Minimum Cut """
#Function to find the Minimum Cut of a Graph
def find_minimum_cut(g):
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
			

    #Return Graph Datastructure without self-loops
	return g

 



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

























