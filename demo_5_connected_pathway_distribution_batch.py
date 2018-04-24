from os import listdir
import math
import matplotlib.pyplot as plt
import networkx as nx

from capstone import Dataclean as dc
from capstone import Subnetwork as sub
from capstone import Parameter as para
from capstone import Stats as st


def main():

    # map loaded
    data_biogrid = dc.biogrid_readfile('BIOGRID-ORGANISM-Homo_sapiens-3.4.156.tab.txt')
    nx_map = dc.to_networkx_map(data_biogrid, 'OFFICIAL_SYMBOL_A', 'OFFICIAL_SYMBOL_B')

    # sample loaded
    sample_current = dc.input_nodelist('test.txt')
    sample_current = dc.nodes_in_map(nx_map, sample_current)

    # pathway loaded
    pathways = {}
    for filename in listdir('pathways'):
        key = filename[0:-4]
        path = 'pathways/' + filename
        pathways[key] = dc.input_nodelist(path)
        pathways[key] = dc.nodes_in_map(nx_map, pathways[key])

    print("You have", len(pathways), 'pathways have been loaded.')

    for k, v in pathways.items():
        print(k)

    for k, v in pathways.items():
        print("Now we're analyzing the pathway: {}".format(k))
        st.connected_pathway_distribution(v, nx_map, num_of_sub=100000, pathway_name=k)


if __name__ == "__main__":
    main()


