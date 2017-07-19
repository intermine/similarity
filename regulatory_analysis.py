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

#Function to create centrality features for the Nodes in the Regulatory Network
def 



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
	




regulatory_analysis()
