""" InterMine @ Open Genome Informatics : 
      -> Implementation of Algorithm depicted in the paper by Blondel : http://epubs.siam.org/doi/pdf/10.1137/S0036144502415960  """
      
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

#Implementation of the Blondel Mathematical Formula for Similarity Matrix Formulation
def perform_blondel_similarity(g):
	#Number of Nodes in the network
	n = len(g.nodes())

	if n < 2:
		print "Error! At least two nodes need to be present in the regulatory network"

		#Exit the function
		return

	#Initial Randomised Matrix
	matrix = np.random.rand(n,n)

	#Adjacency Sparse Matrix
	adjacency_matrix = nx.adjacency_matrix(g)

	#Conversion of Sparse to Dense Matrix
	ad_matrix = adjacency_matrix.todense()

	#Run until convergence to get the similarity matrix
	for i in range(0,100):
		#Formula : Xk+1 = A. Xk . A^T  +  A^T . Xk. A
		matrix = np.matmul(np.matmul(ad_matrix,matrix),ad_matrix.transpose()) + np.matmul(np.matmul(ad_matrix.transpose(),matrix),ad_matrix)


	return matrix


#Function to get the top two similar genes for each gene
def get_top_similar_genes(similarity_matrix,g):
	#Genes in the regulatory network
	genes = g.nodes()

	#Similar Gene List
	similar_gene = []

	#Obtain the top two genes for each gene from the similarity matrix
	for row in similarity_matrix:
		#Conversion to list
		gene_relations = row.tolist()[0]

		#Sort the List according to values
		sorted_relations = sorted(gene_relations,reverse=True)

		#First Similar Match
		first = sorted_relations[0]

		#Second Similar Match
		second = sorted_relations[1]

		#Extraction of indexes for each of the gene
		index1 = gene_relations.index(first)
		index2 = gene_relations.index(second)

		#Add to main similar_gene list
		similar_gene.append([genes[index1],genes[index2]])


	#Similar Dict
	similar_dict = {}

	for i in range(0,len(genes)):
		similar_dict[genes[i]] = similar_gene[i]


	return similar_dict



#Main Function
def main():
	#Connection to Neo4j
	graph = Graph("http://localhost:7474/db/data/cypher",password="rimo")

	#Get a list of Genes and their corresponding length
	genes, lengths = get_genes(graph)

	#Get the regulatory network
	regulatory_network = get_regulatory_networks(graph,genes)

	#Initializing instance for Directed NetworkX Graph
	g = nx.DiGraph()

	#Add the edges into the Graph
	for edge in regulatory_network:
		g.add_edge(edge[0],edge[1])	

	#Call to the function to obtain the similarity matrix
	similarity_matrix = perform_blondel_similarity(g)
    
    #Call to the function to obtain the topmost similar genes for each gene
	top_similar_genes = get_top_similar_genes(similarity_matrix,g)

	

	



main()