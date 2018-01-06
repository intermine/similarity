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

gene_data = [literal_eval(line.strip()) for line in file if line.strip()]

G = nx.Graph() # Initialize the graph

# Adding all the nodes in the Graph and creating edges
for item in gene_data:
	G.add_edge(item[0][0], item[0][1], weight=item[1]) # The data is in the form: ((gene1, gene2), weight)

def get_graph_features(graph):
	# Finding the degree and parity of each node. The degree will help us color the nodes when using D3.js
	print "Calculating the degree and parity for each node"
	for ix, deg in dict(graph.degree()).items():
		graph.node[ix]['degree'] = deg
		graph.node[ix]['parity'] = (1-deg%2)
	
	# Finding Katz centrality of each node. We can use this to decide the width of the links between nodes in the graph 
	print "Generating Katz centrality for each node"
	for ix, katz in nx.katz_centrality(graph, max_iter=max_iterations).items():
		graph.node[ix]['katz'] = katz

	return graph

G = get_graph_features(G)
data = json_graph.node_link_data(G)

print "Writing to 'Visualisations/top_edges_grap.json' file"
with open('Visualisations/top_edges_graph.json', 'w') as file:
	json.dump(data, file, indent=4)


Handler = SimpleHTTPServer.SimpleHTTPRequestHandler

httpd = SocketServer.TCPServer(("", PORT), Handler)

print "serving at port", PORT
httpd.serve_forever()