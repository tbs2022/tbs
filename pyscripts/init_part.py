import numpy as np
import networkx as nx
# from matplotlib import pyplot as plt
import scipy
import igraph as ig
import gurobipy as gp
from gurobipy import GRB

dstdir = "../tmp/"
def runMQP(G):
    print(len(G.nodes),len(G.edges))
    NWeight = [G.nodes[i]['NodeWeight'] for i in G.nodes]
    # EWeight = [G.edges[i]['EdgeWeight'] for i in G.edges]
    nowg = G
    node_list = list(nowg.nodes)
    N = len(node_list)
    W = nx.adj_matrix(nowg,weight='EdgeWeight')

    Wary = W.toarray()
    linesum = np.sum(Wary,0)
    print(linesum)
    maxline = linesum.max() / 2

    num_machines = 12

    #Model
    model = gp.Model('Graph Partition')

    #variables
    map_nodes = model.addMVar((num_machines,N),lb=0.0, ub=1.0,vtype=GRB.BINARY, name="partition")
    max_band = model.addMVar(1,name="max_band")

    #constraint0: each node on one machine only
    certainty = model.addConstrs((map_nodes[:,i].sum() == 1
                                  for i in range(N)), name='certainty')
    #constraint1: CPU etc.
    cap_cpu = 1000
    machine_capacity = model.addConstrs((map_nodes[i,:] @ np.array(NWeight)  <= cap_cpu
                                        for i in range(num_machines)), name='machine_capacity')
    #constraint2: Edge
    from itertools import combinations
    for i,j in combinations(range(num_machines),2):
        psd = map_nodes[i,:] @ W @ map_nodes[j,:] + maxline*sum(map_nodes[k, :] @ map_nodes[k, :] for k in [i,j])
        edge = model.addConstr(psd - max_band <= maxline*(map_nodes[i,:].sum()+map_nodes[j].sum()))

    model.setObjective(max_band,GRB.MINIMIZE)

    model.setParam('Threads', 12)
    model.setParam('NoRelHeurTime', 5)
    model.setParam('TimeLimit', 10)
    # model.setParam('LogFile', logdir+fn+td+"gurobi.log")

    # model.setParam('SolFiles', "/home/data/data12/zhy/tbr/tmp/sols/tmpsol")
    # model.setParam('MIPFocus', 3)

    model.optimize()

    partition = {}
    for j in range(N):
        for i in range(num_machines):
            isin = float(map_nodes[i,j].x)
            if isin > 0.5:
                assert(j not in partition)
                partition[j] = i
    listtowrite = [str(i[1]) + '\n' for i in sorted(partition.items())]
    g = open(dstdir + "initial", "w")
    g.writelines(['#' + str(int(max_band.x)) + '\n'])
    g.writelines(listtowrite)
    g.close()

def readin():
    coars = open(dstdir + "coarse","r")
    alllines = [i.rstrip() for i in coars.readlines()]
    numN = int(alllines[0].split()[0])
    numM = int(alllines[0].split()[1])
    cut = int(alllines[0].split()[2])
    print(numN,numM)
    G = nx.Graph()
    for i in range(1,len(alllines)):
#         G.add_node(i,NodeWeight=1)
        G.add_node(i, NodeWeight=int(alllines[i].split()[0]))

    for i in range(1,len(alllines)):
        thisline = alllines[i].split()
        for j in range(1,len(thisline),2):
            G.add_edge(i,int(thisline[j]), EdgeWeight=int(thisline[j+1]))
#             print(i,thisline[j],int(thisline[j+1]))
#         for j in range(0,len(thisline)):
#             G.add_edge(i,int(thisline[j]),EdgeWeight=1)
    assert(list(G.nodes) == list(range(1,numN+1)))
    runMQP(G)


readin()
