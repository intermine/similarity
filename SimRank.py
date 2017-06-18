""" 
  InterMine @ Open Genome Informatics : Similarity Project
   -> Implementation of the SimRank Algorithm to create a Similarity Matrix for the Gene Regulatory Network
   -> The Similarity Matrix measure will be combined with doc_cluster measure to Rank Genes, in a similar way as to how web pages are ranked        """



#Libraries
from __future__ import division
import json
import pandas as pd
from py2neo import Node, Relationship, Graph
import re
import math
import string
import networkx as nx
from networkx.algorithms.connectivity import minimum_st_edge_cut,minimum_edge_cut
import matplotlib.pyplot as plt
from matplotlib import pylab
import sys
import numpy as np
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.cluster import KMeans,AgglomerativeClustering, MiniBatchKMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.preprocessing import normalize
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from categorical_cluster import hierarchical_mixed
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble,
             discriminant_analysis, random_projection)


#Function to get the Gene Regulatory Network for the corresponding InterMine Model
def get_regulatory_network():
	#Connection to Neo4j
	graph = Graph("http://localhost:7474/db/data/cypher",password="rimo")

	#Extracting Source / Destination
	query = """ MATCH (n:Gene)-[:REGULATES]->(m:Gene) RETURN n.primaryIdentifier,m.primaryIdentifier """
	regulations = graph.data(query)
	
	#Conversion into correct format
	graph = []
	for edge in regulations:		
		graph.append((edge['n.primaryIdentifier'],edge['m.primaryIdentifier']))

	return graph

#Test Function for properties of the Graph
def properties(regulatory_graph):
	#Self Loops
	self_loops = map((lambda graph: graph[0]==graph[1]),regulatory_graph)
	self_loops = sum(self_loops)

    #Returns number of Self-loops
	return self_loops

#Function to populate the Similarity Matrix for the first time
def get_similarity_matrix(di_graph):
	#Shape of Matrix => NXN (N=> Number of Nodes)
	dimension = len(di_graph.nodes())

	#Initializing the numpy matrix with zeros
	matrix = np.zeros((dimension,dimension))

	#Initial Condition : S(a,b) = 1 , for a=b 
	for i in range(0,dimension):
		matrix[i][i] = 1


	return matrix

#Function to Implement the Similarity Score for each pair of incoming edges
def compute_scores(graph,old_matrix,current_matrix,decay,incoming_a,incoming_b,i,j):
	
	#If either of incoming_a or incoming_b is empty -> Set score = 0
	if (not incoming_a) or (not incoming_b):
		current_matrix[i][j] = 0
		return current_matrix

	#Number of Incoming Edges for the first node 
	I_a = len(incoming_a)
	I_b = len(incoming_b)

	#Both of the nodes have at least one incoming edge
	total_score = 0
	for a in incoming_a:
		for b in incoming_b:
			i_a = a[0]
			i_b = b[0]
			first_index = graph.nodes().index(i_a)
			second_index = graph.nodes().index(i_b)
			total_score += old_matrix[first_index][second_index]

	score = (decay / (I_a * I_b)) * total_score

	current_matrix[i][j] = score

	return current_matrix





#Function to perform iterative computation of similarity scores
def calculate_similarity_scores(graph,matrix,iteration,decay):
	#Initial Matrix
	current_matrix = matrix
    
    #Each Iteration will produce a Similarity Matrix which will be used in the next Iteration
	for k in range(0,iteration):
		old_matrix = current_matrix.copy()
		#Calculating S(a,b) = (C/|I(a)| * |I(b)| ) * (For each(i,j) : Sim (I(i),I(j)))
		for i in range(0,len(graph.nodes())):
			for j in range(0,len(graph.nodes())):
				#Node Corresponding to i
				a = graph.nodes()[i]
				#Node Corresponding to j
				b = graph.nodes()[j]

				incoming_a =  graph.in_edges(a)
				incoming_b =  graph.in_edges(b)

				new_matrix = compute_scores(graph,old_matrix,current_matrix,decay,incoming_a,incoming_b,i,j)

				current_matrix = new_matrix


	return current_matrix



#Function to compute the Sim-Rank Matrix
def compute_sim_rank(graph):
	#Conversion to NetworkX Graph
	di_graph = nx.DiGraph()

	#Adding edges
	di_graph.add_edges_from(graph)

	#Node List
	nodes = di_graph.nodes()

	#Conversion into adjacency matrix
	adjacency = nx.to_numpy_matrix(di_graph)

	#Get Initial Similarity Matrix
	similarity_matrix = get_similarity_matrix(di_graph)

	#Final Similarity Matrix
	final_matrix = calculate_similarity_scores(di_graph,similarity_matrix,5,0.5)

	return final_matrix
	





def main():
	#Gene Regulatory Network (Directed Graph)
	regulatory_graph = get_regulatory_network()
	
	#Computing the SimRank Matrix
	similarity_matrix = compute_sim_rank(regulatory_graph)

	#Print
	#print similarity_matrix

	




main()