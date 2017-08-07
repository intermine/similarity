""" 
  InterMine @ Open Genome Informatics : Similarity Project
   -> Clustering the sets by treating each set as a document and clustering them using Document Clustering Methods    
   -> This is a relatively simple model considering no requirement of taking sequence taking into account,
      however the texts are relatively short in nature where tf-idf does not work the best, as for most of the cases : tf*idf ~ idf        
   -> The problem is equivalent to Short Text Clustering                                                                                """


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

from features import create_features, get_genes

#Function to compute the Clusters corresponding the Genetic Sets
def compute_clusters(features):
	#Conversion into numpy array
	features =  np.array(features)

	#As the matrix is sparse with lots of zeros - SVD works better than PCA as it does not centre the data
	svd = TruncatedSVD(n_components=100,random_state=42)

	#Obtaining Features after dimensionality reduction
	transformed_features = svd.fit_transform(features)

	#Normalize the Transformed Feature into Unit Norm so that K-means behave as spherical K-means
	normalized_features = normalize(transformed_features)

	#Range of number of clusters
	cluster_range = 37

	#Silhouette scores
	silhouettes = []

	#Silhouette Analysis
	for cluster_number in range(2,cluster_range):

		#Applying K-means Algorithm - Initialization
		k_means = KMeans(n_clusters=cluster_number, init='k-means++', n_init=10)

		#Fit and get index of each cluster
		cluster_labels = k_means.fit_predict(normalized_features)

		#Average of all the scores of the samples
		silhouette_average = silhouette_score(normalized_features,cluster_labels)

		#print silhouette_average

		#Add to main list
		silhouettes.append(silhouette_average)

	#Obtain max score
	score_max = max(silhouettes)

	#Find index for deciding the number of clusters
	position = silhouettes.index(score_max)

	#Final Number of clusters
	number_of_clusters = position + 2

	#Now apply K-means with the number of clusters determined by Silhouette Analysis
	k_means_final = KMeans(n_clusters=number_of_clusters,init='k-means++',n_init=10)

	#Fitting the Feature Set into Mini Batch K-means
	k_means_final.fit(normalized_features)

	#Cluster Index
	final_index = k_means_final.labels_

	return final_index


#Term frequency for a Genetic Term : Log Normalization
def compute_term_frequency(document,word):
	term_count =  document.count(word)

	if term_count==0:
		return 0 

	return (1 + math.log(term_count))

#Function to compute the inverse document frequency for each document
def compute_idf(gene_documents):
	#Dictionary for storing idf for each word
	idf = {}
	#Vocabulary for the Entire Genetic Information
	vocabulary = list(set([token for document in gene_documents for token in document]))

	#For each token computing idf
	for token in vocabulary:
		#Returns a list consisting of documents where the token is present
		Df = map((lambda document: token in document),gene_documents)

		#Storing token along with Idf in dictionary
		idf[token] = math.log(len(gene_documents) / sum(Df))

	
	return idf

#Function to convert each Gene document into a vector
def compute_tfidf(gene_documents):
	#Inverse Document Frequency for each word
	idf = compute_idf(gene_documents)

	#Final Vector
	vector = []

	for document in gene_documents:
		temp_vector = []
		for term in idf.keys():
			#Tf
			term_frequency = compute_term_frequency(document,term)
			#Idf
			inverse_document_frequency = idf[term]
			#Tf-Idf Score
			tf_idf_score = term_frequency * inverse_document_frequency
			temp_vector.append(tf_idf_score)

		vector.append(temp_vector)

	return vector



#Function to create an equivalent document containing information for each gene
def create_gene_documents(feature):
	#Each Sub-list 
	document_set = []

	#Each Position in the list document_set will correspond to information pertaining to a Gene
	for gene in feature:
		temp = []
		for features in gene:
			temp += features

		document_set.append(temp)

	return document_set	

#Function to obtain the Gene Clusters wrt ID's
def get_gene_clusters(cluster_labels,genes):
	#Get Unique Labels
	unique_labels = list(set(cluster_labels))

	clusters = []

	for cluster_id in unique_labels:
		#Indexes for each cluster
		gene_indexes = np.argwhere(cluster_labels == cluster_id)

		#Extraction from list structure
		gene_clusters = map((lambda x: genes[x[0]]),gene_indexes)

		#Add to main cluster
		clusters.append(gene_clusters)


	return clusters



#Main Function for calls
def main_operation():
	#Connection to Neo4j
	graph = Graph("http://localhost:7474/db/data/cypher",password="rimo")

	#Obtain the features for corresponding InterMine Model
	feature_array = create_features()

	#Obtaining the list of Genes
	genes, length_genes = get_genes(graph)

	#Computing singular sets for each gene as a document =: Type : Lists
	gene_documents = create_gene_documents(feature_array)

	#Conversion into Feature Vectors
	tfidf_vectors = compute_tfidf(gene_documents)

	#Obtaining the cluster labels
	cluster_labels = compute_clusters(tfidf_vectors)

	#Obtain clusters in form of Gene ID's
	gene_clusters = get_gene_clusters(cluster_labels,genes)

	print gene_clusters



main_operation()


