
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Analysis of the Beanmine Dataset                                                      """


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


def main_operation():
	#Connection to Database
	neo4j_graph = Graph('http://localhost:7474/db/data/cypher/',password= 'rimo')




main_operation()