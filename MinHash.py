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



""" Algorithm Description : MinHash  -- This is a method of compressing the information and an approximation scheme towards Jaccard coefficient computation
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
		if gene: #Checking if information pertaining to the Gene is present
			gene_modified = [x for x in gene if x is not None]
			ids = map((lambda g: binascii.crc32(g) & 0xffffffff),gene_modified)
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
			#Selection for coefficient B
			b = hash_b[i]
			#Complete List of Hash Values
			hash_values = map((lambda g: (a*g + b)%c),gene)
			#Choose Minimum Hash Value
			temp.append(min(hash_values))

		minhash_signatures.append(temp)


	return minhash_signatures


""" Application of LSH on the MinHash Signatures to identify similar patterns in the data set 
	Parameters =>   b -> Number of Bands,   r -> Number of Rows in a band,  signature_matrix -> Contains signature for each gene
	                components -> Number of elements in the Gene Signature                                                         """

#Perform LSH
def LSH(b,r,signature_matrix,components):
	# b * r = number of components in a signature
	
	#Initializing LSH matrix - to be filled after hashing bands
	lsh_matrix = []
    
    #Calculating Hash values for each band -> Total Number of hash values for a gene = Number of bands
	for gene in signature_matrix:
		#For storing hash value of each band
		temp = []
		for band in range(0,b):
			#Extracting required number of elements for a band
			temp = gene[band*r:band*r + r]
			#Obtaining a Hash Value
			hash_value = hash(tuple(temp))
			#Append to temporary list
			temp.append(hash_value)

		lsh_matrix.append(temp)

	return lsh_matrix


#Function to churn out the candidate genes which needs to be compared
def candidate_genes(lsh_matrix):
	#Transpose the Matrix 
	transposed_lsh_matrix = np.array(lsh_matrix).transpose()

	#List storing the Candidate Genes from each Band
	candidates = [] 

	#Extracting the candidate genes from each band
	for band in transposed_lsh_matrix:
		#Extract the unique hash values from the band
		hash_list = list(set(band))

		#Candidates for the particular band
		band_candidate = [np.where(band == value) for value in hash_list]
        
        #Add to main candidate list
		candidates.append(band_candidate)


	return candidates


""" The following method taken : N * (N-1) / 2 time, which is not good for scaling  -- Replaced with LSH """
#Function to construct a Similarity Matrix - By comparing all pairs - This will be changed by Locality Sensitive Hashing for faster processing
def get_similarity_matrix(signatures,n,components):
	#Generate Empty Matrix
	similarity_matrix = np.zeros((n,n))

	#Loop through and get common signature values
	for i in range(0,n):
		for j in range(i+1,n):
			#Similarity Metric
			similarity_score = len(set(signatures[i]).intersection(set(signatures[j]))) / len(set(signatures[i]).union(set(signatures[j])))
			similarity_matrix[i][j] = similarity_score


	return similarity_matrix

#Function to get the Similar Genes using information about candidate genes
def get_similar_genes(candidate_gene,genes):
	#Extract Candidate Genes from the Band Hash Values where repetition is observed
	#for band in candidate_gene:
		#print len(band[0][0])
		#print band[0][0]
	return 1




#Main Function for calls
def main(components):
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
	signatures = generate_signatures(shingles,components)

	#similarity_matrix = get_similarity_matrix(signatures,len(genes),components)

	#return similarity_matrix
	
	b = 10 
	r = 4

	#Obtain the matrix formed due to Locality Sensitive Hashing
	lsh_matrix = LSH(b,r,signatures,components)
    
    #Candidate Genes for Close Inspection
	candidate_gene = candidate_genes(lsh_matrix)

	#Use the information regarding candidate genes to obtain similarity scores
	final_similarity = get_similar_genes(candidate_gene, genes)




	




#Compute Similarity Matrix
matrix = main(40)


