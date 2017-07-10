"""  MCL Algorithm   """



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

from features import create_features, get_genes

#Function to find holes in the feature set for the genes
def MCL(graph,inflation,e):
	#Nodes of the Graph
	nodes = graph.nodes()

	#Convert into adjacency matrix
	adjacency = nx.to_numpy_matrix(graph)

	#Transpose the matrix
	adjacency_matrix = adjacency.transpose()

	print adjacency_matrix
	




def main():
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

	inflation = 2
	e = 2

	#Perform MCL Algorithm
	MCL(graph,inflation,e)



	




main()