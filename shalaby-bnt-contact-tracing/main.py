# Explore data from http://www.sociopatterns.org/datasets/infectious-sociopatterns/

import networkx 
import summigraph
path='data/cumulative_2009-04-28.gml';
#G=networkx.read_gml(path);
G=igraph.load(path)
