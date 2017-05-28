
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Application of Cycle Detection Algorithm to understand the dependencies of the Genes  
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


sys.setrecursionlimit(2000)
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



cycle_count = 0

#Function to Implement Depth for Search for finding Cycle
def detect_cycle(start,graph,pred,temp,cycles):
	global cycle_count

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
				
				#print temp
				cycle_count = cycle_count + 1	
				print cycle_count

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
			



	#print graph.edges()




#End of function


#Test Cases for Testing the Algorithm
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





#Loading gene_interactions JSON file into a variable 
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
	if i==300:
		break
	i = i + 1



#print len(graph.nodes())
#print len(graph.edges())

test_graph = nx.Graph()
test_graph = test_cases()
#Creation of appropriate data structure for DFS -- Initially mark all nodes
graph_nodes = {}
for item in graph.nodes():
	graph_nodes[item] = 'w'


#print graph_nodes

#Create Edge List
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



cycle_detection(edge_list)



save_graph(graph,"intermine.pdf")

print cycle_count/2
