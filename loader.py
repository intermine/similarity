""" 
  InterMine @ Open Genome Informatics : Similarity Project
  Loader Script for Neo4j  """



#Libraries 
import json
import pandas as pd
from py2neo import Node, Relationship, Graph
import re
import string


#Loading gene_interactions JSON file into a variable 
with open('JSON rows/gene_interactions.json') as json_data:
	d = json.load(json_data)


#Source List
source = []
#Target List 
target = []
#Storing the results of interactions
result_list = d["results"]

for edge in result_list:
	source.append(edge[0])
	target.append(edge[2])


graph = Graph('http://localhost:7474/db/data/cypher/')

for edge in result_list:
	query = '''MERGE(w0:Gene{gene:{gene1}}) MERGE(w1:Gene{gene:{gene2}}) MERGE (w0)-[:INTERACTS]->(w1)'''
	#graph.run(query,gene1 = edge[0],gene2 = edge[2])



#print len(source)
#print len(target)











