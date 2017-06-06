
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Application of Cycle Detection Algorithm to understand the dependencies of the Genes 
   -> Application of Graph Analysis methods using Neo4j 
   -> Application of Graph Centrality Measures to find Similarity Amongst Nodes
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
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd


sys.setrecursionlimit(10000)
#For Storage of Cycles
cycles = []

#Function for Drawing the Graph
def save_graph(graph,file_name):
	#initialze Figure
	plt.figure(num=None, figsize=(20, 20), dpi=80)
	plt.axis('off')
	fig = plt.figure(1)
	pos = nx.spring_layout(graph)
	nx.draw_networkx_nodes(graph,pos)
	nx.draw_networkx_edges(graph,pos)
	nx.draw_networkx_labels(graph,pos)

	cut = 1.00
	xmax = cut * max(xx for xx, yy in pos.values())
	ymax = cut * max(yy for xx, yy in pos.values())
	plt.xlim(0, xmax)
	plt.ylim(0, ymax)

	plt.savefig(file_name,bbox_inches="tight")
	pylab.close()
	del fig


temp = []
#Initialization
cycle_count = 0

#Function to Implement Depth for Search for finding Cycle
def detect_cycle(start,graph,pred,temp,cycles):
	#For tracking cycle count
	global cycle_count
	cycle_count = cycle_count + 1
	temp.append(start)  
	
	#Color the Visiting Node
	graph_nodes[start]='b'
	#print graph_nodes

	visiting_list = []
	for edge in graph:
		if edge[0] == start:
			visiting_list.append(edge[1])
	#print visiting_list


	for vertex in visiting_list:
		if graph_nodes[vertex]=='b':
			#print "Cycle Found"
			if pred!=vertex:
				
				cycle_count = cycle_count + 1   
				#print cycle_count
				temp = []
			
		else:
			detect_cycle(vertex,graph,start,temp,cycles)
		#temp = []

				
#Function for Detecting Cycles in the Given Network -- Implementation of DFS
def cycle_detection(graph):
	for item in graph_nodes:
		#This means this node is not visited
		if graph_nodes[item] == 'w':
			#Color the node and mark it as visited
			detect_cycle(item,graph,0,[],[])
			


#Test Cases for Testing the Algorithm : Cycle Detection
def test_cases():
	#Temporary Test Graph
	test_graph = nx.Graph() 

	#Test Edges -- Presence of Cycle
	test_graph.add_edge(1,2)
	test_graph.add_edge(2,3)
	test_graph.add_edge(1,3)
	test_graph.add_edge(3,4)
	test_graph.add_edge(3,5)
	test_graph.add_edge(4,5)
	test_graph.add_edge(7,8)
	test_graph.add_edge(7,9)
	test_graph.add_edge(8,9)


	return test_graph

#Function to perform functions via Neo4j operations
def graph_analytics(graph):

	#For finding Triangular cycles in the Graph -- For tracking short range interactions
	triangular_cycle = graph.data("match (a)-[:INTERACTS]->(b)-[:INTERACTS]->(c)-[:INTERACTS]->(a) return distinct a,b,c")
	
	#For Finding Self-Loops in the Graph
	loop = graph.data("match (n)-[r]->(n) return n")

	print triangular_cycle




#Function to get a Path between given a Pair of Nodes
def path_node(graph,node1,node2):
	#Neo4j Query for finding paths  
	query = ''' match p=(n1)-[:INTERACTS*]-(n2) where n1.gene = {gene1} and n2.gene = {gene2} return p'''
	#List of Paths
	paths = graph.data(query,gene1=node1,gene2=node2)


	return paths

#Function to plot the centrality 3-Dimensional Data after PCA
def plot_3D(dataset):
	#Elements along X-axis
	x = [np.take(ele,0) for ele in dataset]
	#Elements along Y-axis
	y = [np.take(ele,1) for ele in dataset]
	#Elements along Z-axis
	z = [np.take(ele,2) for ele in dataset]
	#Conversion into Numpy Array
	x = np.array(x)
	y = np.array(y)
	z = np.array(z)

	fig = plt.figure()
	#Creating a subplot
	ax =fig.add_subplot(111,projection='3d')
	#Scatter Graph with parameters
	ax.scatter(x,y,z,c='r',marker='o')
	#Setting Labels for each dimension
	ax.set_xlabel('X label')
	ax.set_xlabel('Y label')
	ax.set_xlabel('Z label')

	#Displaying the plot
	plt.show()


#Function to perform Silhouette Analysis to determine the optimum amount of clusters
def silhouette_analysis(dataset):
	#Range of number of clusters
	cluster_range = [2,3,4,5,6,7,8,9,10]

	#Silhouette Scores
	silhouette_scores = []

	for cluster in cluster_range:
		#Initialize clusterer with number of clusters as cluster
		clusterer = KMeans(n_clusters=cluster, random_state=10)
		#Fit the dataset and get the index of each cluster
		cluster_labels = clusterer.fit_predict(dataset)
		#Average of the scores of all samples(points)
		silhouette_avg = silhouette_score(dataset, cluster_labels)
		silhouette_scores.append(silhouette_avg)
		print silhouette_avg

	#Obtaining the maximum silhouette score and the number of clusters corresponding to it
	max_score = max(silhouette_scores)
	index = silhouette_scores.index(max_score)

	#Final number of clusters for K-means to run
	final_cluster = cluster_range[index]

	#Cluster corresponding to the optimum cluster number
	clusterer = KMeans(n_clusters=final_cluster,random_state=10)

	#Final labels for the cluster
	labels = clusterer.fit_predict(dataset)

	#print final_cluster

	return labels


""" 
    Description => _function network_centralization
        := Computes Features based on connectivity in Graph and information pertaining to Genes(Pathway,Ontology,Protein Domains)
        := Structural Features based on Centrality Measures
        := Nodal Features are one-hot encoded [Ontology, Domains, Pathway]       """



#Function to find out the most important nodes in the network using Connectivity Measures
def network_centralization(graph,protein_domain,gene_ontology,unique_protein_id,unique_ontologies):
	#Degree Centrality -- Fraction of Node the node is connected to
	centrality_degree = nx.degree_centrality(graph)

	#Closeness Centrality -- Reciprocal of the sum of the short paths from the node to all the other nodes
	centrality_closeness = nx.closeness_centrality(graph)

	#Betweenness Centrality -- Fraction of all pair shortest paths that pass through the node
	centrality_betweenness = nx.betweenness_centrality(graph)

	#Communicability Centrality -- Sum of closed walks of all length starting and ending at node n (Currently commented ->Needs high memory to run)
	#centrality_communicability = nx.communicability_centrality(graph) 

	#Katz Centrality - Number of paths that pass through the node : Variation of closeness centrality
	""" alpha : Attenuation factor  """
	#centrality_katz = nx.katz_centrality(graph,alpha=0.1,beta=1.0,max_iter=20000)

	#PageRank -- Returns the ranking of the nodes based on the structure of links
	""" alpha : Dampening Factor for PageRank """
	page_rank = nx.pagerank(graph,alpha=0.8)


	#Computing a list of centralities for each node
	centralities = {}
	feature_list = []


	#Creation of Dictionary with all centrality measures for each node
	for node in graph.nodes():
		temp = []
		temp.append(centrality_degree[node])
		temp.append(centrality_closeness[node])
		temp.append(centrality_betweenness[node])
		#temp.append(centrality_communicability[node])
		temp.append(page_rank[node])

		#Adding information corresponding to Number of Protein Domains
		try:
			number_of_domains =  len(protein_domain[node])
			temp.append(number_of_domains)
		except:
			temp.append(0)

		#Adding information corresponding to Number of Ontologies
		try:
			ontology_terms = len(gene_ontology[node])
			temp.append(ontology_terms)
		except:
			temp.append(0)

		
		""" Treatment of Gene Ontology Information and Protein Domain Information as categorical variables """
	
		#Adding Protein Domain features 
		try:
			domains = protein_domain[node]
			domain_vector = [0] * len(unique_protein_id)
			for domain in domains:
				position = unique_protein_id.index(domain)
				domain_vector[position] = 1

			temp = temp + domain_vector

		except:
			domain_vector = [0] * len(unique_protein_id)
			temp = temp + domain_vector

		#Adding Gene Ontology Features
		try:
			ontology = gene_ontology[node]
			ontology_vector = [0] * len(unique_ontologies)
			for onto in ontology:
				position = unique_ontologies.index(onto)
				ontology_vector[position] = 1

			temp = temp + ontology_vector

		except:
			ontology_vector = [0] * len(unique_ontologies)
			temp = temp + ontology_vector


		feature_list.append(temp)
		centralities[node] = temp

  
	#Initializing the numpy array
	feature_list = np.array(feature_list)


	pca = PCA(n_components = 3)

	new_list = pca.fit(feature_list)

	#Get the components from transforming the original data
	matrix = pca.transform(feature_list)

	#Plotting the 3D Data
	#plot_3D(matrix)

	#Silhouette Analysis
	node_labels = silhouette_analysis(feature_list)

	return node_labels



""" @Main Function -- Responsible for calling functions which do smaller graph operations """

def main_operation():
	#Loading gene interactions JSON file into a variable 
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
		


	test_graph = nx.Graph()
	test_graph = test_cases()

	#Creation of appropriate data structure for DFS -- Initially mark all nodes
	graph_nodes = {}
	for item in graph.nodes():
		graph_nodes[item] = 'w'


	#Create Edge List for Cycle Detection
	edge_list = []
	for item in graph.edges():
		temp = []
		temp.append(item[0])
		temp.append(item[1])
		edge_list.append(temp)
		temp = []
		temp.append(item[1])
		temp.append(item[0])
		edge_list.append(temp)


	#Protein Domain information
	with open('JSON rows/gene_proteindomains.json') as json_data:
		proteins = json.load(json_data)

	#Storing only essential information
	proteins = proteins["results"]

	#Storing Genes and their corresponding protein domains
	protein_domain = {}

	protein_id = []
	
	#Initial Population
	for protein_info in proteins:
		protein_domain[protein_info[0]] = []
		protein_id.append(protein_info[1])

	for protein_info in proteins:
		protein_domain[protein_info[0]].append(protein_info[1])

	#Unique protein ID's
	unique_protein_id = list(set(protein_id))


	#Gene Ontology Information
	with open('JSON rows/gene_goterms.json') as json_data:
		go_terms = json.load(json_data)


	#Storing only essential information
	go_terms = go_terms["results"]

	#Storing genes and their corresponding GO_terms
	gene_ontology = {}

	#Ontology ID's
	ontology_id = []

	#Initial population
	for onto in go_terms:
		gene_ontology[onto[0]] = []
		ontology_id.append(onto[5])

	for onto in go_terms:
		gene_ontology[onto[0]].append(onto[5])


	#Unique List of Ontologies
	unique_ontologies = list(set(ontology_id))


	#cycle_detection(edge_list)

	#save_graph(graph,"intermine.pdf")

	#print cycle_count/2

	#Initializing variable for Neo4j Analytics
	#neo4j_graph = Graph('http://localhost:7474/db/data/cypher/')
	#Calling function for performing graph analytics on Neo4j
	#graph_analytics(neo4j_graph)
	#Calling function for finding path between two nodes
	#path_node(neo4j_graph,graph.nodes()[0],graph.nodes()[4])

	final_clusters = network_centralization(graph,protein_domain,gene_ontology,unique_protein_id,unique_ontologies)

	#Drawing the Graph
	colors = {}
	for node in graph.nodes():
		index = graph.nodes().index(node)
		colors[node] = final_clusters[index]

	values = [colors.get(node,0.25) for node in graph.nodes()]
	nx.draw(graph,cmap=plt.get_cmap('jet'),node_color = values)
	plt.show()

	#Testing Purpose
	#data = np.random.rand(50,4)
	#silhouette_analysis(data)




#Calling main_operation function for detecting cycles
main_operation()
