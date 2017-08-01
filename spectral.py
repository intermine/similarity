
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Application of the Spectral Clustering Methods to understand similar nodes in the Graph
   -> Treatment of the Data Set as Undirected Graph -- Gene Interactions taken into account """


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
import scipy.sparse as sparse

'''
   Spectral Clustering Algorithm : 
	  -> Based on the connectivity between the nodes
	  Steps: 
			1. Compute Adjacency Matrix for the given Graph
			2. Compute the Laplacian Matrix for the Graph (Diagonal Degree Matrix - Adjacency Matrix)
			3. Normalize the matrix
			4. Get the Eigen Values and Eigen Vectors corresponding to the Normalized Laplacian Matrix
			5. Choose the number of clusters (say k), then choose k eigen vectors corresponding to eigen values
			6. Apply K-means on them to cluster the nodes with number of clusters as k

	  @_parameters : 
		 lap => normalized laplacian matrix
		 eigen_number => The number of Eigen vectors to be obtained

	  @_return :
		 A numpy array with cluster labels corresponding to each node 
'''

#Function to compute the eigen vectors and eigen value & Apply K-means
def spectral_clustering(lap,eigen_number):

	#Eigen Values and Eigen Vectors for the Normalized Laplacian Matrix
	w,v = np.linalg.eigh(lap)

	#Convert the Eigen Vector to list 
	eigen_vectors = v.tolist()

	#Extraction of the indexes of the top few eigen values
	indexes = sorted(range(len(w)), key=lambda j: w[j])[-10:]
    
    #Feature vector for the top eigen vectors
	feature_vector =  v[:,indexes]	
	
	#Apply K-means to cluster the nodes
	clustered = KMeans(n_clusters=4,random_state=10)

	#Predict the Labels
	cluster_labels = clustered.fit_predict(v)

	#Returns a numpy array consisting of cluster assignment to each node
	return cluster_labels



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

		if i == 200:
			break

		i += 1

	#Converting into numpy array
	adjacency_matrix = nx.to_numpy_matrix(graph)

	#Checking if the Adjacency Matrix Generated is symmetric
	#print np.allclose(adjacency_matrix,adjacency_matrix.T,atol=1e-8)

	#Normalized Laplacian Matrix
	normalized_laplacian = nx.normalized_laplacian_matrix(graph)

	#Normal Laplacian Matrix
	laplacian = nx.laplacian_matrix(graph)
    
    #Conversion to dense matrix
	dense_laplacian = laplacian.todense()

	#Conversion into Dense Matrix
	normalized_laplacian_dense = normalized_laplacian.todense()

	#The number of eigen vectors to be obtained for K-means clustering
	eigen_number = 2

	#Function Call for Spectral Clustering
	cluster_labels = spectral_clustering(normalized_laplacian_dense,eigen_number)

	#List of Nodes
	nodes = graph.nodes()

	#Cluster consisting of similar nodes
	clusters = {}

	#Initializing the clusters
	for num in range(0,eigen_number):
		clusters[num] = []

	i = 0

	#Populating with the actual node names in the cluster
	for cluster_number in cluster_labels:
		clusters[cluster_number].append(nodes[i])
		i+=1

	
	#Visualization
	cluster_dict = {}
	for node in clusters[0]:
		cluster_dict[node] = 0

	for node in clusters[1]:
		cluster_dict[node] = 1


	#Drawing the Graph
	values = [cluster_dict.get(node,0.25) for node in graph.nodes()]
	#nx.draw(graph,cmap=plt.get_cmap('jet'),node_color = values)
	#plt.show()



	


#Call the main function for Implementing Spectral Clustering 
main_operation()