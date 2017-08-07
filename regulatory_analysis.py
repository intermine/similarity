#Regulatory Graph Analysis

#Libraries
from __future__ import division
import json
import pandas as pd
from py2neo import Node, Relationship, Graph
import re
import math
import string
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import pylab
import sys
import numpy as np
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.cluster import KMeans,AgglomerativeClustering, MiniBatchKMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.preprocessing import normalize
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble,
             discriminant_analysis, random_projection)
from features import get_regulatory_networks, get_genes
from tensor import AutoEncoder
from sklearn.mixture import GaussianMixture
from sklearn.manifold import Isomap

#Function to create centrality features for the Nodes in the Regulatory Network
def create_features(graph):
	#Nodes in the regulatory network
	genes = graph.nodes()

	#In Degree Centrality
	in_degree = nx.in_degree_centrality(graph)

	#Out Degree Centrality
	out_degree = nx.out_degree_centrality(graph)

	#Closeness Centrality
	closeness = nx.closeness_centrality(graph)

	#Betweeness centrality -- Fraction of all-pairs shortest path that passes through a given network
	node_betweenness = nx.betweenness_centrality(graph)

	#Page Rank for node ranking values
	page_rank = nx.pagerank(graph,alpha=0.8)

	#Dictionary for storing features for each gene
	feature_dictionary = {}

	#Populate with Genes and Features
	for gene in genes:
		#Initialization
		feature_dictionary[gene] = []

		#In Degree
		feature_dictionary[gene].append(in_degree[gene])

		#Out Degree
		feature_dictionary[gene].append(out_degree[gene])

		#Closeness
		feature_dictionary[gene].append(closeness[gene])

		#Node Betweenness
		feature_dictionary[gene].append(node_betweenness[gene])

		#Page Rank
		feature_dictionary[gene].append(page_rank[gene])


	return feature_dictionary

#Function to find important edges based on centrality - Edge Betweenness
def important_edges(graph):
	#Edge Betweenness - To find important edges in the network
	edge_betweenness = nx.edge_betweenness_centrality(graph)

	return edge_betweenness


#Function to create visualization
def visualization(features,cluster_labels):
	#Elements along X-axis
	x = [np.take(ele,0) for ele in features]

	#Elements along Y-axis
	y = [np.take(ele,1) for ele in features]

	fig = plt.figure()

	plt.scatter(x,y,c=cluster_labels,label='')
	#plt.scatter(x,y,label='')
	plt.show()

#Function to create new features from trained network
def create_new_features(feature_array,weights,biases):
	#New feature list
	features = []

	#Append to empty Numpy Array
	for feature in feature_array:
		features.append(np.add(np.matmul(feature,weights),biases).tolist())
     
	return np.array(features)


#Non Linear Dimensionality Reduction using Isomaps - Manifold Learning
def iso_map(feature_array):
	#Parameters : Number of neighbours to be considered for the point graphs and number of components
	isomap = Isomap(n_neighbors = 4,n_components = 2)

	#Components in lower dimensional space
	components = isomap.fit_transform(feature_array)

	return components



#PCA 
def principal_component_analysis(feature_array):
	pca = PCA(n_components = 2)
	pca.fit(feature_array)
	#visualization(pca.transform(feature_array),pca.transform(feature_array))

	return pca.transform(feature_array)


#Function to cluster the features
def cluster(features):
	#Range of number of clusters
	cluster_range = [2,3,4,5,6]

	#Silhouette Scores
	silhouette_scores = []

	for cluster in cluster_range:
		#Initialize clusterer with number of clusters as cluster
		clusterer = KMeans(n_clusters=cluster,random_state=10)
		#Fit the dataset and get the index of each cluster
		cluster_labels = clusterer.fit_predict(features)
		#Average of the scores of all samples(points)
		silhouette_avg = silhouette_score(features, cluster_labels)
		silhouette_scores.append(silhouette_avg)
		print silhouette_avg

	#Average Silhouette Score to be set as threshold
	average_silhouette_score = sum(silhouette_scores) / float(len(silhouette_scores))

	#Extraction of indices for cluster values above silhouette threshold
	cluster_above_threshold = [silhouette_scores.index(i)+2 for i in silhouette_scores if i > average_silhouette_score]

	return cluster_above_threshold


#Function to compare the results of two methods
def compare_nonlinear(isomap_labels,autoencoder_labels):
	#Total labels
	label_length = len(isomap_labels)

	#Percentage
	percentage = 0

	#Clusters for isomap
	clusters_isomap = []

	#Clusters due to autoencoders
	clusters_autoencoders = []

	for i in range(0,3):
		clusters_isomap.append(np.where(isomap_labels == i)[0])
		clusters_autoencoders.append(np.where(autoencoder_labels == i)[0])


	#Conversion into appropriate format
	clusters_isomap = [item.tolist() for item in clusters_isomap]
	clusters_autoencoders = [item.tolist() for item in clusters_autoencoders]	
    

	for i in range(0,label_length):
		#Element is i
		for j in range(0,3):
			if i in clusters_isomap[j] and i in clusters_autoencoders[j]:
				percentage += 1

	print percentage / label_length
	



#Function to perform Expectation Maximization
def perform_EM(new_features):
	#Estimation of Model Parameters through Expectation maximization results
	gmm = GaussianMixture(n_components=3, covariance_type='full').fit(new_features)

	cluster_labels = gmm.predict(new_features)

	return cluster_labels


#Function to Extracted the Top Edges from the Network
def print_top_edges(top_edges):
	#Top Sorted Edges
	top_edges_sorted = sorted(top_edges.items(),key = lambda x: -x[1])

	#Save the Results in a File
	f_edges = open('Results/top_edges.txt','w')

	for edge in top_edges_sorted:
		f_edges.write(str(edge))
		f_edges.write("\n\n")


	f_edges.close()




#Function to create features for nodes in the Regulatory Graph
def regulatory_analysis():
	#Connection to Neo4j
	graph = Graph("http://localhost:7474/db/data/cypher",password="rimo")

	#Get a list of Genes and their corresponding length
	genes, lengths = get_genes(graph)

	#Get the regulatory network
	regulatory_network = get_regulatory_networks(graph,genes)

	#Initializing instance for Directed NetworkX Graph
	g = nx.DiGraph() 

	#Conversion into list of tuples
	edge_list = map((lambda x : (x[0],x[1])),regulatory_network)

	#Addition of edges into the directed graph
	g.add_edges_from(edge_list)

	#Creation of features for the nodes
	node_features = create_features(g)

	#Feature List
	feature_list = []

	for gene in node_features:
		feature_list.append(node_features[gene])


	#Conversion into numpy array
	feature_array = np.array(feature_list)

	#Autoencoder for Dimensionality Reduction
	weights, biases = AutoEncoder(feature_array)

	#Create 2D features using trained weights and biases
	new_features = create_new_features(feature_array,weights,biases)

	#Clustering - Generating clusters above threshold
	#cluster_above_threshold = cluster(new_features)

	#Obtain Labels from EM Algorithm for Autoencoder features
	gaussian_labels_autoencoders = perform_EM(new_features)	

	#Visualization
	#visualization(new_features,gaussian_labels)

	##iso_map_features = iso_map(feature_array)
	#visualization(new_features,new_features)
	

	#Obtain Labels from EM Algorithm for Iso map features
	##gaussian_labels_isomap = perform_EM(iso_map_features)

	



	#visualization(iso_map_features,gaussian_labels_isomap)

	##visualization(new_features,gaussian_labels_autoencoders)

	##pca_feature = principal_component_analysis(feature_array)

	##gaussian_labels_pca = perform_EM(pca_feature)

	##compare_nonlinear(gaussian_labels_isomap,gaussian_labels_pca)

	#visualization(pca_feature,gaussian_labels_pca)
	#visualization(,gaussian_labels_isomap)

	#top_edges = important_edges(g)

	#print_top_edges(top_edges)

	









regulatory_analysis()
