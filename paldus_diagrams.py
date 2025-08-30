from params import *
import paldus_classes
from paldus_cas import cas
from paldus_cas import swap_count
from paldus_cas import cas_to_ugg

from paldus_basic import *
from paldus_classes import ugg
from copy import deepcopy
from paldus_classes import arithmetic_string
from paldus_classes import disambiguate
from collections import deque
from paldus_classes import pair_permutations
import math
from itertools import product
from itertools import permutations
from paldus_classes import integrate
from paldus_classes import virtual, occupied, general
from fortran_code import *
from templates import *
import sys
import time
import pickle
import pydot
from pydot import Dot, Edge, Node
# from graphviz import Graph



def simple_diag():
    print('witam')
    graph = pydot.Dot(graph_type='graph')
    # let's add the relationship between the king and vassals

    for i in range(3):
        # we can get right into action by "drawing" edges between the nodes in our graph
        # we do not need to CREATE nodes, but if you want to give them some custom style
        # then I would recomend you to do so… let's cover that later
        # the pydot.Edge() constructor receives two parameters, a source node and a destination
        # node, they are just strings like you can see
        edge = pydot.Edge("king", "lord%d" % i)
        # and we obviosuly need to add the edge to our graph
        graph.add_edge(edge)
        
    # now let us add some vassals
    vassal_num = 0
    for i in range(3):
        # we create new edges, now between our previous lords and the new vassals
        # let us create two vassals for each lord
        for j in range(2):
            edge = pydot.Edge("lord%d" % i, "vassal%d" % vassal_num)
            graph.add_edge(edge)
            vassal_num += 1
            
    # ok, we are set, let's save our graph into a file
    graph.write_png('example1_graph.png')
        
    # and we are done!




    # creating nodes is as simple as creating edges!
    # http://www.graphviz.org/doc/info/attrs.html
    # which in turn is part of the full docs in
    # http://www.graphviz.org/Documentation.php

    # neat, huh? Let us create the rest of the nodes!
    

    
    node_s1 = pydot.Node("S1", style="filled", fillcolor='lightgray', rank="same")
    node_s2 = pydot.Node("S2", style="filled", fillcolor='lightgray', rank="same")
    node_s3 = pydot.Node("S3", style="filled", fillcolor='lightgray')
    node_s4 = pydot.Node("S4", style="filled", fillcolor='lightgray')

    node_t1 = pydot.Node("T1", style="filled", fillcolor="lightblue")
    node_t2 = pydot.Node("T2", style="filled", fillcolor="lightblue")
    node_t3 = pydot.Node("T3", style="filled", fillcolor="lightblue")
    node_t4 = pydot.Node("T4", style="filled", fillcolor="lightblue")


    #ok, now we add the nodes to the graph
    graph = pydot.Dot(rank="same")
    graph.add_node(node_s1)#, rank="same")
    graph.add_node(node_s2)
    graph.add_node(node_s3)
    graph.add_node(node_s4)

    # .add_node(node_t1)
    # S.add_node(node_t2)
    # S.add_node(node_t3)
    # S.add_node(node_t4)

    S = pydot.Subgraph(rank='same')
    S.add_node(node_t1)
    S.add_node(node_t2)
    S.add_node(node_t3)
    S.add_node(node_t4)

    D = pydot.Subgraph(rank='same')
    D.add_node(node_s1)
    D.add_node(node_s2)
    D.add_node(node_s3)
    D.add_node(node_s4)

    # graph.add_node(node_t1)
    # graph.add_node(node_t2)
    # graph.add_node(node_t3)
    # graph.add_node(node_t4)

    graph.add_edge(pydot.Edge(node_s1, node_t1, label='k'))
    graph.add_edge(pydot.Edge(node_t1, node_s2, label='c'))
    graph.add_edge(pydot.Edge(node_t2, node_s3, label='e'))
    graph.add_edge(pydot.Edge(node_s4, node_t4, label='m'))
    graph.add_edge(pydot.Edge(node_t4, node_s4, label='d'))
    graph.add_edge(pydot.Edge(node_s3, node_t3, label='l'))

    graph.add_edge(pydot.Edge(node_t2, node_t3, style="invis"))
    graph.add_edge(pydot.Edge(node_s2, node_s3, style="invis"))
    graph.add_edge(pydot.Edge(node_s1, node_s2, arrowhead='none'))
    graph.add_edge(pydot.Edge(node_s3, node_s4, arrowhead='none'))
    graph.add_edge(pydot.Edge(node_t1, node_t2, arrowhead='none'))
    graph.add_edge(pydot.Edge(node_t3, node_t4, arrowhead='none'))
    graph.add_subgraph(D)
    graph.add_subgraph(S)
    # graph.add_edge(pydot.Edge(node_a, node_b))

    graph.write_png('example3_graph.png')


    node1 = pydot.Node(1)
    node2 = pydot.Node(2)
    node3 = pydot.Node(3)
    node4 = pydot.Node(4)
    
    P = pydot.Dot()
    P.add_edge(pydot.Edge(node1,node2))
    P.add_edge(pydot.Edge(node2,node3))
    P.add_edge(pydot.Edge(node1,node4))
    
    S = pydot.Subgraph(rank='same')
    S.add_node(node3)
    S.add_node(node4)
    P.add_subgraph(S)
    P.write_png('foo.png')
