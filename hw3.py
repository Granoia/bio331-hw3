import lab3
import random
import math


#bio331 homework 3


def RWR(adj_ls, start_node, q, t):
    counts = {}
    for node in adj_ls:
        counts['node'] = 0
    
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
        counts[e] = math.log(counts[e])
        



def main():
    #some stuff
    node_ls, edge_ls, edge_dict = lab3.readData('EGFR1-reachable.txt')
    adj_ls = lab3.make_adj_ls(node_ls, edge_ls)
    



if __name__ == '__main__':
    main()
