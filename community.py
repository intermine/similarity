
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Implementation of the Louvain Algorithm for Community Detection -- Using Fast Heuristics 
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





def apply_louvain(nodes,edges):
	#Initialization for sum of weights of all links
	m = 0
	#Sum of all weights incident on a given node
	incident_weight = {}
	for node in nodes:
		incident_weight[node] = 0

	for edge in edges:
		m += edge[2]
		incident_weight[edge[0]] += edge[2]
		#incident_weight[edge[1]] += edge[2]

	#Halved : Undirected Nature of Graph
	m = m/2

	#Initial Step : Assigning Different Community to all the Different Nodes
	best_partition = [[node] for node in nodes]
	partition = first_phase(nodes,edges)




def first_phase(nodes,edges):
	#Calling partition function for initially Assigning Communities
	best_partition = initial_partition(nodes,edges)



def second_phase(nodes,edges):
	return 1


def initial_partition(nodes,edges):
	#Initially assigning communities to nodes
	partition = [[node] for node in nodes]




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

	#Function Call to apply the Louvain Algorithm
	apply_louvain(nodes,edges)


	



main()