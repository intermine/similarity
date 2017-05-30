
"""  InterMine @ Open Genome Informatics : Similarity Project
   -> Application of the Spectral Clustering Methods to understand similar nodes in the Graph
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


