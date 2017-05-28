
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

#Function to Implement Depth for Search for finding Cycle
def detect_cycle(start,graph,current,temp,cycles):
	#Checking the adjacent edge and recursively calling the function
	for edge in graph.edges():
		if start == edge[0]:
			end = edge[1]
			for item in graph_nodes:
				if item[0] == end:
					color = item[1]
					#Already Visited Node -- Cycle Found
					if color == 'b':
						print "Cycle Found"
						cycles.append(temp)
						print temp
						#print current

						cycles.append(1)

					else:
						item[1] ='b'
						current = end
						temp.append(end)
						detect_cycle(end,graph,current,temp,cycles)
					break

		temp = []


				

			

#Function for Detecting Cycles in the Given Network -- Implementation of DFS
def cycle_detection(graph):
	for item in graph_nodes:
		#This means this node is not visited
		if item[1] == 'w':
			#Color the node and mark it as visited
			item[1] = 'b'
			current = item[0]
			detect_cycle(item[0],graph,current,[],[])




#End of function


#Test Cases for Testing the Algorithm
def test_cases():
	#Temporary Test Graph
	test_graph = nx.Graph()

	#Test Nodes
	test_graph.add_node(1)
	test_graph.add_node(2)
	test_graph.add_node(3)
	test_graph.add_node(4)
	test_graph.add_node(5)
	test_graph.add_node(6)
	test_graph.add_node(7)
	test_graph.add_node(8)
	test_graph.add_node(9)

	#Test Edges
	test_graph.add_edge(1,2)
	test_graph.add_edge(2,3)
	test_graph.add_edge(2,4)
	test_graph.add_edge(4,5)
	test_graph.add_edge(4,6)
	test_graph.add_edge(5,6)
	test_graph.add_edge(7,8)
	test_graph.add_edge(7,9)
	test_graph.add_edge(8,9)
	test_graph.add_edge(3,5)


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

#Extracting the edge relationships
for edge in graph_info:
	source.append(edge[0])
	target.append(edge[2])

	#Adding the edge in NetworkX
	graph.add_edge(edge[0],edge[2])



#print len(graph.nodes())
#print len(graph.edges())

test_graph = nx.Graph()
test_graph = test_cases()
#Creation of appropriate data structure for DFS -- Initially mark all nodes
graph_nodes = []
for item in test_graph.nodes():
	temp = []
	temp.append(item)
	temp.append('w')
	graph_nodes.append(temp)

cycle_detection(test_graph)


#save_graph(graph,"intermine.pdf")

