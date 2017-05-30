
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Implementation of the Louvain Algorithm for Community Detection -- Using Fast Heuristics 
   -> Treatment of the Data Set as Undirected Graph -- Just Interactions taken into account """



#Libraries 
import community
import json
import pandas as pd
from py2neo import Node, Relationship, Graph
import re
import string
import networkx as nx
from networkx.algorithms.connectivity import minimum_st_edge_cut,minimum_edge_cut
import matplotlib.pyplot as plt



#Declarations for Functions
self_loops = {}
edges_of_nodes = {}
communities = []
actual_partition = []

def apply_louvain(nodes,edges):
	#Global variable for self-loops
	global self_loops
	
	#Initialization for sum of weights of all links
	m = 0
	#Sum of all weights incident on a given node
	incident_weight = {}
	for node in nodes:
		incident_weight[node] = 0

	#Initially remove any self-edges if present in the graph -- Initially the algorithm would not be able to handle self loops
	for edge in edges:
		if edge[0] == edge[1]:
			edges.remove(edge)

	#Edges corresponding to a particular Node
	global edges_of_nodes
	

	#Communities
	global communities 
	communities = [node for node in nodes]
	
	global actual_partition
	actual_partition = []

	#Storing the edges corresponding to each Node
	for edge in edges:
		m += edge[2]
		incident_weight[edge[0]] += edge[2]
		#incident_weight[edge[1]] += edge[2]

		#For each node add the edges corresponding to it
		if edge[0] not in edges_of_nodes:
			edges_of_nodes[edge[0]] = [edge]
		else:
			edges_of_nodes[edge[0]].append(edge)


	

	#Halved : Undirected Nature of Graph
	m = m/2

	#End of Initialization

	network = (nodes,edges)

	#Initial Step : Assigning Different Community to all the Different Nodes
	best_partition = [[node] for node in nodes]

	partition = first_phase(network)



	



#Function to perform first phase of louvain algorithm
def first_phase(network):
	#Calling partition function for initially Assigning Communities
	best_partition = initial_partition(network)




#Function to get the Neighbours
def neighbours(node):
	neighbour_node = []
	for item in edges_of_nodes[node]:
		#Self-Loop
		if item[0]==item[1]:
			continue
		else:
			neighbour_node.append(item[1])

	return neighbour_node




#Function for Initializing Partition
def initial_partition(network):
	#Initially assigning communities to nodes
	nodes = network[0]
	edges = network[1]
	partition = [[node] for node in nodes]

	#Initialize self_loops
	for node in nodes:
		self_loops[node] = 0 
	
	#Adding information about self-loops over here
	for edge in edges:
		#Presence of Self-Loop
		if edge[0]==edge[1]:
			self_loops[edge[0]] += edge[2]
			self_loops[edge[1]] += edge[2]

	#Return a list of nodes according to partition
	return partition





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

	G = nx.erdos_renyi_graph(30, 0.05)

	#first compute the best partition
	partition = community.best_partition(G)



	



main()