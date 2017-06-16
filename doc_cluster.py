""" 
  InterMine @ Open Genome Informatics : Similarity Project
   -> Main script for calling other functions in the module
   -> Finds Similarity amongst sets                                        """


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

from features import create_features
from categorical_cluster import hierarchical_mixed


def main_operation():
	#Obtain the features for FlyMine Model
	feature_array = create_features()
	









main_operation()


