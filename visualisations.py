from networkx.readwrite import json_graph
import networkx as nx
import codecs
from ast import literal_eval
import json
import SimpleHTTPServer
import SocketServer

max_iterations = 10000 # This is so that the katz_centrality function can converge for every node.
PORT = 8000 # SimpleHTTPServer PORT

file = codecs.open('Results/top_edges.txt', encoding='utf-8')

# Load gene data and weights
print "loading gene data..."
gene_data = [literal_eval(line.strip()) for line in file if line.strip()]

# This function is used to normalize the data.
def normalize(data_list):
	min_val = min(data_list)
	max_val = max(data_list)
	return [ float(item-min_val)/float(max_val-min_val) for item in data_list ]

# The available weights are normalized so that they can be used to show the strength of the link between two nodes
weights = normalize([item[1] for item in gene_data])

G = nx.Graph() # Initialize the graph

# Adding all the nodes in the Graph and creating edges
for i in range(len(gene_data)):
	# The data is in the form: ((gene1, gene2), original_weight)
	G.add_edge(gene_data[i][0][0], gene_data[i][0][1], weight=weights[i]) 

def get_graph_features(graph):
	# Finding the degree of each node. The degree will help us color the nodes when using D3.js
	print "Calculating the degree and parity for each node..."
	for ix, deg in dict(graph.degree()).items():
		graph.node[ix]['degree'] = deg
		
	return graph


G = get_graph_features(G)

data = json_graph.node_link_data(G, {'source': 'source', 'target': 'target', 'id': 'id'})

print "Writing to 'Visualisations/top_edges_grap.json' file..."
with open('Visualisations/top_edges_graph.json', 'w') as file:
	json.dump(data, file, indent=4)

Handler = SimpleHTTPServer.SimpleHTTPRequestHandler

httpd = SocketServer.TCPServer(("", PORT), Handler)

print "serving files at: localhost:"+ str(PORT) + "/Visualisations/top_edges.html"

httpd.serve_forever()
