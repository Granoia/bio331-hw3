#!/usr/bin/python

## utility functions for writing JSON files.
## BIO331
## Anna Ritz
import json

def test():
	"""
	Test function to make sure that the module is successfully imported.
	"""
	print('json_utils properly imported!')
	return

def write_json(data,jsonfile):
	"""
	Writes the data object as a JSON file.

	:param data: dictionary from make_json_data() function.
	:param jsonfile: string -- name of JSON file.
	"""
	print('\nWriting JSON for graph to outfile %s' % (jsonfile))
	json.dump(data,open(jsonfile,'w'),indent=4)
	return

def make_json_data(nodes,edges,node_attributes=None,edge_attributes=None,title="",description="",tags=[],labels=True):
	"""
	Creates a dictionary that contains the following entries::

		{
			"graph": {
				"nodes": [],
				"edges": []
			},
			"metadata": {
				"title": "",
				"description": "",
		 		"tags": []
			}
		}

	See http://ec2-52-41-252-78.us-west-2.compute.amazonaws.com/help/programmers for more information about the JSON structure.

	:param nodes: list or set of nodes
	:param edges: list or set of edges
	:param node_attributes: dictionary of node attributes. Optional.  The key is a node name and the value is a dictionary of "attribute"/"value" pairs.Dictionary has the form::

		{"node_name" : {"attribute1":"val1","attribute2":"val2",...}}

	:param edge_attributes: dictionary of edge attributes. Optional.  The key is the source node and the value is a dictionary where the key is the target node and the value is a dictionary of "attribute"/"value" pairs.  Dictionary has the form::

		{"source" : {"target" : {"attribute1":"val1","attribute2":"val2",...}}}

	:param title: string -- title of graph. Optional.
	:param description: string -- description of graph. Optional.
	:param tags: list -- list of tag names. Optional.
	:param labels: boolean -- If True, write node IDs as the text on the node.  Default is True.  Set labels=False to not write the names on the nodes.
	:returns: dictionary formatted for JSON.
	"""
	data = {}

	## add metadata.
	data['metadata'] = {'title':title,'description':description,'tags':tags}

	## add graph dictionary with empty lists for nodes and edges.
	data['graph'] = {'nodes':[],'edges':[]}

	## add each node to the data dictionary.
	for node_name in nodes:
		## create node_element dictionary.  id is required.
		node_element = {'id':node_name}

		## if labels=True, add the content attribute to write the node name.
		if labels:
			node_element['content'] = node_name

		## if other attributes are specified, add them in bulk with an "update" function.
		if node_attributes != None and node_name in node_attributes:
			node_element.update(node_attributes[node_name])
		
		## final node for JSON is written as "data" as the key and the 
		## node_element dictionary as the value.
		node_wrapper = {'data':node_element}

		## append node_wrapper to the list of nodes.
		data['graph']['nodes'].append(node_wrapper)

	## add each edge to the data dictionary.
	for source,target in edges:
		## create edge_element dictionary.  source and target are required.
		edge_element = {'source':source,'target':target}

		## if other attributes are specified, add them in bulk with an "update" function.
		if edge_attributes != None and source in edge_attributes and target in edge_attributes[source]:
			edge_element.update(edge_attributes[source][target])

		## final edge for JSON is written as "data" as the key and the
		## edge_element dictionary as the value.
		edge_wrapper = {'data':edge_element}

		## append edge_wrapper to the list of edges.
		data['graph']['edges'].append(edge_wrapper)

	return data