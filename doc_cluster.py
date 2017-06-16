""" 
  InterMine @ Open Genome Informatics : Similarity Project
   -> Clustering the sets by treating each set as a document and clustering them using Document Clustering Methods    
   -> This is a relatively simple model considering no requirement of taking sequence taking into account,
      however the texts are relatively short in nature where tf-idf does not scale well          
   -> The problem is equivalent to Short Text Clustering                                             """


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
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans,AgglomerativeClustering 
from sklearn.metrics import silhouette_samples, silhouette_score
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from categorical_cluster import hierarchical_mixed
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble,
             discriminant_analysis, random_projection)

from features import create_features, get_genes

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



#Main Function for calls
def main_operation():
	#Connection to Neo4j
	graph = Graph("http://localhost:7474/db/data/cypher",password="rimo")

	#Obtain the features for FlyMine Model
	feature_array = create_features()

	#Obtaining the list of Genes
	genes, length_genes = get_genes(graph)

	#Computing singular sets for each gene as a document
	gene_documents = create_gene_documents(feature_array)

	print gene_documents

	



	









main_operation()


