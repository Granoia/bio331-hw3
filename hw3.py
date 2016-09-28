#########
######### BIO 331 HOMEWORK 3
######### NICK FRANZESE
#########



import random
import math
import json_utils
import graphspace_utils


############################################
#Data reading functions (there's only one)##
############################################

def readData(filename):
    """
    reads in the data, returns a list of nodes, a list of edges, and a dictionary whose keys are edges and whose values are Interaction Types.
    """
    node_set = set()
    edge_ls = []
    edge_type_dict = {}
    with open(filename, 'r') as f:
        for line in f:
            data = line.strip('\n').split('\t')
            node1 = data[0]
            node2 = data[1]
            node_set.add(node1)
            node_set.add(node2)
            if node1 != node2:
                edge = (data[0], data[1])
                edge_type_dict[edge] = data[2]
                edge_ls.append(edge)

    return set_to_list(node_set), edge_ls, edge_type_dict


def set_to_list(s):
    """
    helper function for readData() and probably some other things
    takes a set, returns a list with the contents of the set
    """
    ls = []
    for item in s:
        ls.append(item)
    return ls
    
    

    

############################################
#Graph processing functions#################
############################################


def make_adj_ls(nodes, edges):
    """
    make an adjacency list in the form of a dictionary whose keys are nodes and whose values are dictionaries
    the subdictionaries have keys 'out' for a list of outgoing edges, 'visited' for visited status, and 'distance' for distance from start node (the latter two are for BFS)
    """
    d = {}
    for n in nodes:
        d[n] = {}
        d[n]['out'] = []        #list of outgoing edges
        d[n]['visited'] = False #visited status for BFS
        d[n]['distance'] = None #distance from start node for BFS
    
    for e in edges:
        node1 = e[0]
        node2 = e[1]
        if e[0] != e[1]:
            d[node1]['out'].append(node2)

    return d


#accession functions for the adj_ls
def get_neighbors(adj_ls, node):
    return adj_ls[node]['out']

def get_visited(adj_ls, node):
    return adj_ls[node]['visited']

def visit(adj_ls, node):
    adj_ls[node]['visited'] = True

def devisit(adj_ls, node):
    adj_ls[node]['visited'] = False

def reset_visits(adj_ls):
    for k in adj_ls:
        devisit(adj_ls, k)

def get_dist(adj_ls, node):
    return adj_ls[node]['distance']
        
def set_dist(adj_ls, node, d):
    adj_ls[node]['distance'] = d
        

    
def BFS_distances(adj_ls, start_node):
    """
    alters the adj_ls so that distances and visited status reflect the results of the BFS from start_node
    """
    Q = queue()
    Q.enqueue(start_node)
    visit(adj_ls, start_node)
    set_dist(adj_ls,start_node,0)
    while Q.length() != 0:
        w = Q.dequeue()
        w_dist = get_dist(adj_ls,w)
        for n in get_neighbors(adj_ls, w):
            if not get_visited(adj_ls, n):
                visit(adj_ls, n)
                Q.enqueue(n)
                set_dist(adj_ls, n, w_dist + 1)

    return

#queue for implementation of BFS_distances()
class queue:
    def __init__(self):
        self.q = []

    def enqueue(self,new):
        self.q.append(new)

    def dequeue(self):
        if len(self.q) == 0:
            return None
        else:
            ret = self.q.pop(0)
            return ret

    def length(self):
        return len(self.q)




def RWR(adj_ls, start_node, q, t):
    """
    random walk with restarts
    returns a normalized dictionary of the frequency with which a node was walked on
    """
    counts = {}
    for node in adj_ls:
        counts[node] = 0
    
    current = start_node
    for step in range(t):
        counts[current] += 1
        r = random.uniform(0,1)
        if r <= q:
            if len(adj_ls[current]['out']) == 0:
                current = start_node
            else:
                current = random.choice(adj_ls[current]['out'])
        else:
            current = start_node

    for e in counts:
        if counts[e] != 0:
            counts[e] = math.log(counts[e])
            
    new = count_normalizer(counts)
        
    return new


    
    
    
    
    
##########################################################################################
#Functions that process the results of the graph processing data (mostly normalization)###
########################################################################################## 
    
    

def count_normalizer(counts):
    """
    helper function for RWR
    takes a log transformed counts dictionary, returns a normalized version of that dictionary
    """
    new_dict = {}
    max_count = 0
    for c in counts:
        if counts[c] > max_count:
            max_count = counts[c]
    
    min_count = max_count
    
    for c in counts:
        if counts[c] < min_count:
            min_count = counts[c]
            
    for c in counts:
        norm = 1 - ((counts[c] - min_count) / (max_count - min_count))
        new_dict[c] = norm
    
    return new_dict
    
    


    
    
def BFS_d_normalizer(adj_ls):
    """
    takes an adj_ls THAT HAS ALREADY BEEN MODIFIED BY BFS (that's important)
    returns a dictionary whose keys are nodes and whose values are normalized distances
    """
    new_dict = {}
    max_d = get_max_dist(adj_ls)
    
    for n in adj_ls:
        new_dict[n] = (adj_ls[n]['distance']/float(max_d))
    
    return new_dict


def get_max_dist(adj_ls):
    #helper function for BFS_d_normalizer()
    max_dist = 0
    for node in adj_ls:
        if max_dist < get_dist(adj_ls,node):
            max_dist = get_dist(adj_ls,node)
    return max_dist


    
def norm_of_diffs(norm_dist, norm_counts):
    """
    makes a dictionary that normalizes the positive and negative values when you take the difference between the normalized differences and the normalized counts.
    so the node with -1 ends up being the node with the most relative difference favoring RWR and the node with 1 ends up being the node with the most relative difference favoring shortest path.
    """
    new_dict = {}
    neg_ls = []
    pos_ls = []
    for n in norm_dist:
        diff = norm_dist[n] - norm_counts[n]
        if diff < 0:
            neg_ls.append(n)
        else:
            pos_ls.append(n)
        new_dict[n] = diff

    neg_max = 0
    for n in neg_ls:
        if abs(new_dict[n]) > abs(neg_max):
            neg_max = new_dict[n]

    pos_max = 0
    for n in pos_ls:
        if new_dict[n] > pos_max:
            pos_max = new_dict[n]

    for n in neg_ls:
        new_dict[n] = new_dict[n]/neg_max * -1

    for n in pos_ls:
        new_dict[n] = new_dict[n]/pos_max

    return new_dict #now normalized on the negative side and positive side

    
    
    
    
    
    
###########################################################################
#graphspace attribute functions############################################
###########################################################################

def getEdgeAttributes(edges,edge_type_dict=None):
    attrs = {}
    for e in edges:
        source = e[0]
        target = e[1]
        if source not in attrs:
            attrs[source] = {}
        attrs[source][target] = {}
        attrs[source][target]['width'] = 2
        attrs[source][target]['target_arrow_shape'] = 'triangle'
        attrs[source][target]['target_arrow_fill'] = 'filled'
        attrs[source][target]['target_arrow_color'] = 'black'
        if edge_type_dict != None:
            attrs[source][target]['popup'] = '<b>Edge Type:</b> ' + edge_type_dict[e]
    return attrs


def getRWRNodeAttributes(nodes, counts):
    """
    gets the node attributes for the RWR graph
    """
    attrs = {}
    for n in nodes:
        attrs[n] = {}
        attrs[n]['id'] = n
        attrs[n]['content'] = n
        attrs[n]['background_color'] = counts_to_color(counts, n)
    return attrs

    
def getCompNodeAttributes(nodes, norm_dist, norm_counts):
    """
    gets the node attributes for the RWR / shortest path comparison graph
    """
    attrs = {}
    norm_diff_dict = norm_of_diffs(norm_dist,norm_counts)
    for n in nodes:
        attrs[n] = {}
        attrs[n]['id'] = n
        attrs[n]['content'] = n
        attrs[n]['background_color'] = diff_to_color(n, norm_diff_dict)
    return attrs


def counts_to_color(counts, node):
    """
    takes a normalized dictionary of counts, converts it to a shade of green
    """
    if 0 <= counts[node] < .0001: #since the normalized dictionary contains floats, this deals with float 0
        return '#{:02x}{:02x}{:02x}'.format(255,255,255)
    
    
    g = int(counts[node] * 255)
    return '#{:02x}{:02x}{:02x}'.format(0,g,0)
    
def diff_to_color(node, norm_diff_dict):
    """
    finds the difference between the RWR and shortest path distribution for a given node and converts it into a color.
    if it's comparatively more traveled from the RWR, the node will be majenta. If the path is comparatively more distant, it will be green. If it's equally traveled and distant, the node will be white
    """
    norm_diff = norm_diff_dict[node]
    
    if -0.0001 <= norm_diff <= 0.0001:
        return '#{:02x}{:02x}{:02x}'.format(255,255,255)
        
    elif norm_diff < 0:
        modifier = int(norm_diff * -1 * 255)
        return '#{:02x}{:02x}{:02x}'.format(255,255-modifier,255)
    
    elif norm_diff > 0:
        modifier = int(norm_diff * 255)
        return '#{:02x}{:02x}{:02x}'.format(255-modifier,255,255-modifier)

def main():
    node_ls, edge_ls, edge_dict = readData('EGFR1-reachable.txt')
    adj_ls = make_adj_ls(node_ls, edge_ls)
    BFS_distances(adj_ls, 'EGF')
    norm_dist = BFS_d_normalizer(adj_ls)
    norm_counts = RWR(adj_ls, 'EGF', 0.9, 1000000)
    edge_attrs = getEdgeAttributes(edge_ls, edge_dict)
    rwr_attrs = getRWRNodeAttributes(node_ls, norm_counts)
    comp_attrs = getCompNodeAttributes(node_ls, norm_dist, norm_counts)

    data = json_utils.make_json_data(node_ls,edge_ls,rwr_attrs,edge_attrs,'Homework 3 - Random Walk with Restarts','Desc.',['Tag'])
    json_utils.write_json(data,'hw3-1.json')
    graphspace_utils.postGraph('hw3-1','hw3-1.json','franzni@reed.edu','bio331')
    
    
    data2 = json_utils.make_json_data(node_ls, edge_ls,comp_attrs,edge_attrs,'Homework 3 - Difference Between RWR and Shortest Path','Desc.',['Tag'])
    json_utils.write_json(data2,'hw3-2.json')
    graphspace_utils.postGraph('hw3-2','hw3-2.json','franzni@reed.edu','bio331')


#I can't really discern any difference when I change the number of timesteps, but that might be because my normalization method on the differences graph is weird?
#I have a feeling that there should be less and less difference between the two graphs as the number of timesteps increases, but since I normalized that makes it difficult to tell.
#
#Changing q definitely makes a difference though. The nodes that are further away gain a greater amount of difference as q goes down (since it's harder for the random walk to reach the outer nodes when it's more likely
#to restart.


    
if __name__ == '__main__':
    main()
