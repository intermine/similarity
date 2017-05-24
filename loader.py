""" 
  InterMine @ Open Genome Informatics : Similarity Project
   -> Data Cleaning and Extraction
   -> Loads essential Data into Neo4j   """



#Libraries 
import json
import pandas as pd
from py2neo import Node, Relationship, Graph
import re
import string


#Loading gene_interactions JSON file into a variable 
with open('JSON rows/gene_interactions.json') as json_data:
	interactions = json.load(json_data)


#Source List
source = []
#Target List 
target = []
#Storing the results of interactions
result_list = interactions["results"]
#Storing the relevant information - Source , Target, Mode of Interaction (Interaction Detail Type + Interaction Detail Name)
relevant_edge_info = []



#Extracting the node & edge relationships and storing it temporarily - Separately -  for future use
for edge in result_list:
	source.append(edge[0])
	target.append(edge[2])

#Graph Instance for connecting to Neo4j
graph = Graph('http://localhost:7474/db/data/cypher/')

#Loading each edge into Neo4j
for edge in result_list:
	query = '''MERGE(w0:Gene{gene:{gene1}}) MERGE(w1:Gene{gene:{gene2}}) MERGE (w0)-[:INTERACTS]->(w1)'''
	#graph.run(query,gene1 = edge[0],gene2 = edge[2])


""" Extraction and Storage of Protein Domain Information  """

#Loading gene_proteindomains information
with open('JSON rows/gene_proteindomains.json') as json_data:
	domains = json.load(json_data)

#Storage of the essential part 
domains = domains["results"]

#Dictionary for storage of domain information corresponding to each gene
protein_domain_info = {}

#Initialization - With respect to Gene Primary Identifier
for domain_info in domains:
	protein_domain_info[domain_info[0]] = []

#Populating the domains corresponding to each gene
for domain_info in domains:
	temp_domain = []
	#First Element : Domain Identification Number
	temp_domain.append(domain_info[1])
	#Second Element : Domain Name
	temp_domain.append(domain_info[2])
	protein_domain_info[domain_info[0]].append(temp_domain)


""" Extraction and Storage of Pathway Information """

#Loading Pathway Information
with open('JSON rows/gene_pathways.json') as json_data:
	pathways = json.load(json_data)

#Storage of essential part
pathways = pathways["results"]

#Dictionary for storage of pathway information corresponding to each gene
pathway_info = {}

#Initialization - With respect to Gene Primary Identifier
for path in pathways:
	pathway_info[path[2]] = []


#Populating the pathways into the list
for path in pathways:
	temp_list = []
	#First Element : Pathway Identifier
	temp_list.append(path[4])
	#Second Element : Pathway Name
	temp_list.append(path[5])
	#Storage format  : 'gene_ID' : [[pathway_id,pathway_name],[pathway_id,pathway_name],[pathway_id,pathway_name]....]
	pathway_info[path[2]].append(temp_list)


""" Extraction and Storage of Gene Ontology Information """

#Loading Gene Ontology Information
with open('JSON rows/gene_goterms.json') as json_data:
	ontology = json.load(json_data)


#Storage of Essential Part
ontology = ontology["results"]

ontology_info = {}

#Initialization
for onto_info in ontology:
	ontology_info[onto_info[0]] = []


#Population of Ontology Information for each gene
for onto_info in ontology:
	temp_onto = []
	#First Element : OntologyTerm Identifier
	temp_onto.append(onto_info[5])
	#Second Element : OntologyTerm Name
	temp_onto.append(onto_info[6])
	ontology_info[onto_info[0]].append(temp_onto)


""" Storage of Edge Information with regard to Gene Interactions """ 

#Extracting relevant edge information
for edge in result_list:
	temp_edge_list = []
	#First Element : Node1 [Source / Target]
	temp_edge_list.append(edge[0])
	#Second Element : Node2  [Source / Target]
	temp_edge_list.append(edge[2])
	#Third Element : Interaction type
	temp_edge_list.append(edge[4])
	#Fourth Element : Interaction Name
	temp_edge_list.append(edge[5])
	#Final Append
	relevant_edge_info.append(temp_edge_list)


""" Attaching relationships with Edges for entry into Neo4j """

#Creation of appropriate properties for loading the relationships into Neo4j -- Currently treating edges as undirected
for edge in relevant_edge_info:
	#Source 
	source = edge[0]
	#Target
	target = edge[1]

	#Protein Domains for the Source Node
	try:
		source_domain = protein_domain_info[source]
	except:
		source_domain = []

	#Protein Domains for the Target Node
	try:
		target_domain = protein_domain_info[target]
	except:
		target_domain = []


	#Pathways involved with Source Node 
	try:
		source_pathway = pathway_info[source]
	except:
		source_pathway = []

	#Pathways involved with Target Node
	try:
		target_pathway = pathway_info[target]
	except:
		target_pathway = []



	#Ontologies Involved with Source Node
	try:
		source_ontology = ontology_info[source]
	except:
		source_ontology = []

	#Ontologies Involved with Target Node
	try:
		target_ontology = ontology_info[target]
	except:
		target_ontology = 



	








	
















