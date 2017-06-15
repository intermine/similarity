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
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans,AgglomerativeClustering 
from sklearn.metrics import silhouette_samples, silhouette_score
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from categorical_cluster import hierarchical_mixed
from matplotlib import offsetbox
from sklearn import (manifold, datasets, decomposition, ensemble,
             discriminant_analysis, random_projection)


#Function to get the Entire Gene Regulatory Network for Community Detection
def get_regulatory_networks(graph,genes):
	#Extracting Source / Destination
	query = """ MATCH (n:Gene)-[:REGULATES]->(m:Gene) RETURN n.primaryIdentifier,m.primaryIdentifier """
	regulations = graph.data(query)
	
	#Conversion into correct format
	graph = []
	for edge in regulations:
		temp = []
		temp.append(edge['n.primaryIdentifier'])
		temp.append(edge['m.primaryIdentifier'])
		graph.append(temp)

	return graph
	

#Function to extract the Genes that regulate a given Gene
def regulatory_networks_regulated(graph,genes):
	#For each gene, the regulated genes are extracted
	regulated_genes = {}

	for gene in genes:
		regulated_query = """ MATCH (n:Gene{primaryIdentifier:{gene_id}})<-[:REGULATES]-(r) RETURN r.primaryIdentifier  """
		regulates = graph.data(regulated_query,gene_id=gene)
		regulated_genes[gene] = []
		for reg in regulates:
			regulated_genes[gene].append(reg['r.primaryIdentifier'])

	return regulated_genes


#Function to extract the Genes regulated by a given Gene
def regulatory_networks_regulates(graph,genes):
	#For each gene, the regulated genes are extracted
	regulated_genes = {}

	for gene in genes:
		regulated_query = """ MATCH (n:Gene{primaryIdentifier:{gene_id}})-[:REGULATES]->(r) RETURN r.primaryIdentifier  """
		regulates = graph.data(regulated_query,gene_id=gene)
		regulated_genes[gene] = []
		for reg in regulates:
			regulated_genes[gene].append(reg['r.primaryIdentifier'])

	return regulated_genes


#Function to get the chromosome in which the gene is present
def chromosome(graph,genes):
	#For each gene, the chromosome identifier
	chromosomes = {}

	for gene in genes:
		chromosome_query = """ MATCH (n:Gene{primaryIdentifier:{gene_id}})-[:LOCATED_ON]->(chr) RETURN chr.primaryIdentifier """
		gene_chromosome = graph.data(chromosome_query,gene_id=gene)
		chromosomes[gene] = []
		for ch in gene_chromosome:
			chromosomes[gene].append(ch['chr.primaryIdentifier'])

	return chromosomes


#Function the extract the pathways associated with each gene
def pathway(graph,genes):
	#For each gene , the Pathway Identifier
	pathways = {}

	for gene in genes:
		pathway_query = """ MATCH (n:Gene{primaryIdentifier:{gene_id}})<-[:CONTAINS_GENE]-(path) RETURN path.identifier """
		gene_pathway = graph.data(pathway_query,gene_id=gene)
		pathways[gene] = []
		for path in gene_pathway:
			pathways[gene].append(path['path.identifier'])

	return pathways

#Function to extract the phenotypes associated with each gene
def phenotype(graph,genes):
	#For each gene, the Phenotype ID is associated
	phenotypes = {}

	for gene in genes:
		phenotype_query = """ MATCH (n:Gene{primaryIdentifier:{gene_id}})-[:HAS_PHENOTYPE]->(ph) RETURN ph.reagentId """
		gene_phenotype = graph.data(phenotype_query,gene_id=gene)
		phenotypes[gene] = []
		for ph in gene_phenotype:
			phenotypes[gene].append(ph['ph.reagentId'])

	return phenotypes


#Function to extract the tissue in which the gene is expressed 
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


def get_genes(graph):
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

	return genes, lengths


#Function to extract information from relationships and generate a feature array for each Gene
def create_features():
	#Connection to Neo4j
	graph = Graph("http://localhost:7474/db/data/cypher",password="rimo")

	#Get a list of Genes and their corresponding length
	genes, lengths = get_genes(graph)

	""" Extraction of Relationships for feature constructions """

	#Protein Domains
	protein_domains = domains(graph,genes)

	#Diseases Associated with a Gene
	tissues = tissue(graph,genes)

	#Phenotypes Associated with a Gene
	phenotypes = phenotype(graph,genes)

	#Pathways
	pathways = pathway(graph,genes)

	#Chromosomes
	chromosomes = chromosome(graph,genes)

	#Regulates
	regulates = regulatory_networks_regulates(graph,genes)

	#Regulated By
	regulated_by = regulatory_networks_regulated(graph,genes)

	#Regulatory Graph
	regulatory_graph = get_regulatory_networks(graph,genes)

	#Feature Creation
	feature_array = []

	for gene in genes:
		temp = []
		#Protein Domains
		temp.append(protein_domains[gene])
		#Tissues
		temp.append(tissues[gene])
		#Phenotypes
		temp.append(phenotypes[gene])
		#Pathways
		temp.append(pathways[gene])
		#Chromosomes
		temp.append(chromosomes[gene])
		#Regulates
		temp.append(regulates[gene])
		#Regulated By
		temp.append(regulated_by[gene])

		feature_array.append(temp)


	return feature_array




#Function call for Feature Creation
result = create_features()
