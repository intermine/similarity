""" 
  InterMine @ Open Genome Informatics : Similarity Project
   -> Extraction of essential data from Neo4j and generation of feature arrays   """


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

#Function to extract the tissue in which they are expressed 
def tissue(graph,genes):
	#For each gene, the tissues associated with them
	tissues = {}

	for gene in genes:
		tissue_query = """ MATCH (n:Gene{primaryIdentifier:{gene_id}})-[:EXPRESSED_IN]->(tissue) RETURN tissue.name """
		gene_tissue = graph.data(tissue_query,gene_id=gene)
		tissues[gene]=[]
		for tis in gene_tissue:
			tissues[gene].append(tis['tissue.name'])


	return tissues



#Function to extract the Protein Domains for a corresponding Gene
def domains(graph,genes):
	#For each gene, list of protein domains will be stored
	protein_domains = {}

	for gene in genes:
		domain_query = """ MATCH (n:Gene{primaryIdentifier:{gene_id}})<-[:DOMAIN_HAS]-(domain) RETURN domain.name """
		domain = graph.data(domain_query,gene_id=gene)
		protein_domains[gene] = []
		for prot_domain in domain:
			protein_domains[gene].append(prot_domain['domain.name'])


	return protein_domains





#Function to extract information from relationships and generate a feature array for each Gene
def create_features():
	#Establishing Connection
	graph = Graph("http://localhost:7474/db/data/cypher",password="rimo")

	#Query to get the Genes along with their length
	query = "MATCH (n:Gene) RETURN n.primaryIdentifier,n.length LIMIT 10"
	result = graph.data(query)

	#Gene List
	genes = [gene['n.primaryIdentifier'] for gene in result]
	#Length -- Some list fields do not have any value - Hence handle them through exception
	lengths = []
	for gene in result:
		try:
			lengths.append(gene['n.length'])
		except:
			lengths.append(0)

	""" Extraction of Relationships for feature constructions """

	#Protein Domains
	#protein_domains = domains(graph,genes)

	#Diseases Associated with a Gene
	tissues = tissue(graph,genes)


	


	



create_features()