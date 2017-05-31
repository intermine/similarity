
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Application of the Spectral Clustering Methods to understand similar nodes in the Graph
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
from matplotlib import pylab
import sys
import numpy as np 
import collections
from sklearn.cluster import KMeans

#Function to compute the eigen vectors and eigen value & Apply K-means
def spectral_clustering(lap):
	#Eigen Values and Eigen Values for the Normalized Laplacian
	#w, v = np.linalg.eig(lap)
	#s = np.asanyarray(lap)
	

	print lap[0].size



	

	return 1

#Base function for making function calls
def main_operation():
	#Loading gene interactions JSON file into a variable 
	with open('JSON rows/gene_interactions.json') as json_data:
		interactions = json.load(json_data)


	#Information about Graph Connectivity
	graph_info = interactions["results"]

	#Creating NetworkX instance
	graph = nx.Graph()

	i = 0
	#Extracting the edge relationships
	for edge in graph_info:
		#Adding the edge in NetworkX
		graph.add_edge(edge[0],edge[2])

	#Converting into numpy array
	adjacency_matrix = nx.to_numpy_matrix(graph)

	#Conversion into Laplacian matrix
	#laplacian_matrix = nx.laplacian_matrix(graph)

	#Normalized Laplacian Matrix
	normalized_laplacian = nx.normalized_laplacian_matrix(graph)

	#Function Call for Spectral Clustering
	spectral_clustering(normalized_laplacian)





	



main_operation()