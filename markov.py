
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Application of Highly Connected Sub-graph  [HCS] algorithm on the Preliminary Dataset 
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



print len(graph.nodes())
print len(graph.edges())


save_graph(graph,"intermine.pdf")

