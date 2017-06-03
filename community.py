
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Implementation of the Girvan Newman Algorithm for Community Detection 
   -> Girvan Newman : Hierarchial - Divisive in nature : Begin with one group containing all points and then divide successfully
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


#Function to compute the degree of the nodes
def get_degree(graph,adjacency_matrix):
	#To be returned by the function
	degree_dict = {}

	#Node List
	node_list = graph.nodes()

	#Total Number of Nodes
	number_nodes = len(graph.nodes())

	#Computing the sum of each row
	Sum = adjacency_matrix.sum(axis=1)

	for i in range(0,number_nodes):
		degree_dict[node_list[i]] = Sum[i,0]

	#Return a dictionary consisting of nodes and their corresponding degree
	return degree_dict







''' Girvan Newmann Algorithm : 
      Basic overview : 
             Edge Betweenness => Fraction of all pair shortest path that pass through a Given Edge
             Bridge => Considered to be an edge which has a high edge Betweenness
             B/w distinct communities there should be presence of bridges

      Steps: 
            1. Compute the edge Betweenness for the graph
            2. Choose the edge with the highest edge betweenness and remove into
            3. Go to step 1 until no edges are left (If modularity is improving then add to the community structure)


      Measure of checking Community Structure : Modularity
                                                                                                   
                                                                                                   '''	



#Function to extract the edge with maximum edge betweenness and remove it from the graph
def girvan_newman(graph):
	initial_components = nx.number_connected_components(graph)
	temp_comp = initial_components
	print temp_comp
	while temp_comp<=initial_components:
		#Edge Betweenness values -- Returns dictionary of edges with centrality as values
		edge_betweenness = nx.edge_betweenness_centrality(graph)
		#Finding edge with maximum betweenness
		max_edge = max(edge_betweenness.values())
		for k,v in edge_betweenness.iteritems():
			if v == max_edge:
				#Remove the edge from the network
				graph.remove_edge(k[0],k[1])

		#Recalculating the number of connected components -- Will break out of the loop once an extra component is formed due to removal of edge
		temp_comp = nx.number_connected_components(graph)







#Function to call girvan_newman until and unless all edges are consumed -- Checked with respect to modularity at every step
def compute_girvan_newman(graph,degree_node):
	return 1

	


	






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
		graph.add_edge(edge[0],edge[2],weight=1.0)
		if i == 2000:
			break
		i = i + 1

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

	#Adjacency Matrix -- Returns a Scipy Sparse Matrix
	adjacency_matrix = nx.adj_matrix(graph)

	#Number of Nodes
	n = nx.number_of_nodes(graph)

	#Calculating the weighted degree for each node
	degree_node = get_degree(graph,adjacency_matrix)

	#print degree_node

	#m = 0.0
	#for i in range(0,n):
#		for j in range(0,n):
	#		m += adjacency_matrix[i,j]


	#print m


	compute_girvan_newman(graph,degree_node)

	



	



main()