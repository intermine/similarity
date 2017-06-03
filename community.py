
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Implementation of the Girvan Newmann Algorithm for Community Detection -- Using Fast Heuristics 
   -> Treatment of the Data Set as Undirected Graph -- Just Interactions taken into account """



#Libraries 
import json
import pandas as pd
from py2neo import Node, Relationship, Graph
import re
import string
import networkx as nx
from networkx.algorithms.connectivity import minimum_st_edge_cut,minimum_edge_cut
import matplotlib.pyplot as plt







	






def main():
	#Loading gene_interactions JSON file into a variable 
	with open('JSON rows/gene_interactions.json') as json_data:
		interactions = json.load(json_data)


	#Information about Graph Connectivity
	graph_info = interactions["results"]

	source = []
	target = []

	#Creating NetworkX instance
	graph = nx.Graph()

	i = 0
	#Extracting the edge relationships
	for edge in graph_info:
		#temp = []
		source.append(edge[0])
		target.append(edge[2])

		#Adding the edge in NetworkX
		graph.add_edge(edge[0],edge[2])

	#Nodes :[node1,node2,node3.......]
	nodes = []
	for node in graph.nodes():
		nodes.append(node)

	#Edges : [[node1,node2,weight],[],[]....]
	edges = []
	for edge in graph.edges():
		temp = []
		#For undirected conditions
		temp.append(edge[0])
		temp.append(edge[1])
		temp.append(1)
		edges.append(temp)
		temp = []
		temp.append(edge[1])
		temp.append(edge[0])
		temp.append(1)
		edges.append(temp)

	



	



main()