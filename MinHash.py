""" 
  InterMine @ Open Genome Informatics : Similarity Project
   -> Implementation of the MinHash Algorithm to finding similarity amongst Genes treated as sets of information
   -> Implementation of LSH on the signatures generated to find highly similar genes                                     """


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
import binascii
import numpy as np
import pandas as pd
from features import create_features, get_genes
from doc_cluster import create_gene_documents

#Function to Map the Shingles to an ID - using CRC32 Hash
def generate_shingle_id(sets):
	#List containing Shingle ID's
	shingle_ids = []
    
    #Finding the shingle ID's for each set of property in the Gene Set
	for gene in sets:
		ids = map((lambda g: binascii.crc32(g) & 0xffffffff),gene)
		shingle_ids.append(ids)


	return shingle_ids


#Main Function for calls
def main():
	#Connection to Neo4j
	graph = Graph("http://localhost:7474/db/data/cypher",password="rimo")

	#Obtain the features for corresponding InterMine Model
	feature_array = create_features()

	#Obtaining the list of Genes
	genes, length_genes = get_genes(graph)

	#Treating each Gene as a set
	sets = create_gene_documents(feature_array)

	#Singles mapped to ID's
	shingles = generate_shingle_id(sets)

	





main()


