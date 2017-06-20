""" 
  InterMine @ Open Genome Informatics : Similarity Project
   -> Implementation of the MinHash Algorithm to finding similarity amongst Genes treated as sets of information
   -> Implementation of LSH on the signatures generated to find highly similar genes      
   -> This method of Finding has a huge benefit when dealing with millions of sets of Genes to find Similarity as it scales very well          """


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
import random



""" Algorithm Description : MinHash 
       -> Generate Random Hash functions of the form :  (a * x  +  b ) mod c , where  a,b < max(x) and c is the prime number just greater than x
       -> x can take maximum of 2**32 - 1 as in the earlier step for Shingle ID's we have used crc32  
       -> Choose each hash function and generate hash values for each Shingle ID -- then choose the minimum value as signature                    
       -> In this way, if there are 'n' Hash Functions there will be n components of a signature                                                """


#Function to Map the Shingles to an ID - using CRC32 Hash
def generate_shingle_id(sets):
	#List containing Shingle ID's
	shingle_ids = []
    
    #Finding the shingle ID's for each set of property in the Gene Set
	for gene in sets:
		ids = map((lambda g: binascii.crc32(g) & 0xffffffff),gene)
		shingle_ids.append(ids)


	return shingle_ids

#Function to Generate Random values for a & b in the Hash function with the requirement that there is no repetition
def generate_random_hash(max,components):
	#A sampling method cannot be used to generate a list of random number with no repetitions as they don't work in low memory conditions
	#print list(set(np.random.randint(11,size=10)))

	initial_list = np.random.randint(max,size=components)

	while len(set(initial_list)) != len(initial_list):
		#Perform the operation again to get unique random values
		initial_list = np.random.randint(max,size=components)


	return initial_list	


#Function to Generate the MinHash Signatures using Hash Functions
def generate_signatures(shingles,components):
	#Max Shingle_id 
	max_shingle_id = 2**32 - 1

	#Prime Number just greater than max_shingle_id -- Can be extracted from http://compoasso.free.fr/primelistweb/page/prime/liste_online_en.php
	c = 4294967311

	#Generate Random Hash Functions
	hash_a = generate_random_hash(max_shingle_id,components)
	hash_b = generate_random_hash(max_shingle_id,components)
    
    #Storing the Min Hash Signatures in a list
	minhash_signatures = []

	for gene in shingles:
		temp = []
		for i in range(0,components):
			#Selection for coefficient A
			a = hash_a[i]
			#Choose for coefficient B
			b = hash_b[i]
			#Complete Hash values
			hash_values = map((lambda g: (a*g + b)%c),gene)
			#Choose minimum hash value
			temp.append(min(hash_values))

		minhash_signatures.append(temp)


	return minhash_signatures



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

	#Get signature based on the Shingle ID's
	signatures = generate_signatures(shingles,10)

	print signatures





main()


