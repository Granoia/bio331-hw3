import graphspace_utils
import json_utils

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


def handle_complexes(nodes, edges):
    new_node_set = set()
    new_edge_set = set()
    edge_type_dict = {}
    for n in nodes:
        if ',' in n:
            c = n.split(',')
            for item in c:
                new_node_set.add(item)
        else:
            new_node_set.add(n)

    for e in edges:
        if e[0] == e[1]:
            pass
        else:
            node1 = e[0]
            node2 = e[1]
            handle_edge_help(node1,node2,new_edge_set)

    return set_to_list(new_node_set), set_to_list(new_edge_set)


def handle_edge_help(node1,node2,s):
    """
    helper function for handle_complexes that breaks down complexes and adds the appropriate pairwise combinations to the edge set.
    """
    if ',' in node1:
        c1 = node1.split(',')
    else:
        c1 = [node1]

    if ',' in node2:
        c2 = node2.split(',')
    else:
        c2 = [node2]

    for a in c1:
        for b in c2:
            tup = (a,b)
            s.add(tup)



def get_int_types(etd):
    """
    gets the interaction types from the edge type dictionary and returns them
    """
    int_types = set()
    for k in etd:
        int_types.add(etd[k])
    return int_types


def set_to_list(s):
    """
    takes a set, returns a list with the contents of the set
    """
    ls = []
    for item in s:
        ls.append(item)
    return ls


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
            

def get_max_dist(adj_ls):
    #helper function for dist_to_color()
    max_dist = 0
    for node in adj_ls:
        if max_dist < get_dist(adj_ls,node):
            max_dist = get_dist(adj_ls,node)
    return max_dist
    
def dist_to_color(adj_ls, node, max_dist):
    d = get_dist(adj_ls, node)
    if d == 0: #start node will be white, furthest nodes will be completely blue
        return '#{:02x}{:02x}{:02x}'.format(255,255,255)

    rel_dist = float(d)/max_dist
    maxHexValue = 255
    b = int(rel_dist * maxHexValue)
    return '#{:02x}{:02x}{:02x}'.format(0,0,b)
    


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


def getNodeAttributes(nodes, adj_ls):
    attrs = {}
    max_dist = get_max_dist(adj_ls)
    for n in nodes:
        attrs[n] = {}
        attrs[n]['id'] = n
        attrs[n]['content'] = n
        attrs[n]['background_color'] = dist_to_color(adj_ls, n, max_dist)
    return attrs


def main():
    nodes, edges, edge_dict = readData('EGFR1-reachable.txt')
    adj_ls = make_adj_ls(nodes, edges)
    BFS_distances(adj_ls,'EGF')
    
    edgeAttrs = getEdgeAttributes(edges,edge_dict)
    nodeAttrs = getNodeAttributes(nodes, adj_ls)

    #data = json_utils.make_json_data(nodes,edges,nodeAttrs,edgeAttrs,'Lab 3','Desc.',['Tag'])
    #json_utils.write_json(data,'lab3.json')
    #graphspace_utils.postGraph('lab3-2-1','lab3.json','franzni@reed.edu','bio331')


    hc_nodes, hc_edges = handle_complexes(nodes, edges)
    hc_adj_ls = make_adj_ls(hc_nodes,hc_edges)
    BFS_distances(hc_adj_ls,'EGF')
    hc_edgeAttrs = getEdgeAttributes(hc_edges)
    hc_nodeAttrs = getNodeAttributes(hc_nodes, hc_adj_ls)

    

    
    data2 = json_utils.make_json_data(hc_nodes, hc_edges,hc_nodeAttrs,hc_edgeAttrs,'Lab 3: Handle Complexes','Desc.',['Tag'])
    json_utils.write_json(data2,'lab3-hc.json')
    graphspace_utils.postGraph('lab3-4','lab3-hc.json','franzni@reed.edu','bio331')
    
    #print(get_int_types(edge_dict))

    print(len(nodes), len(hc_nodes))
    print(len(edges), len(hc_edges))




if __name__ == '__main__':
    main()
