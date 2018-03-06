# Project: BMI Capstone Project
# Filename: capstone.py
# Time: 201801
# Editor: Hsiu-Ping Lin

import csv
import random
import numpy as np
import math
from matplotlib import pyplot as plt
import pandas as pd
from collections import Counter
import networkx as nx


class Dataclean(object):
    """ This class contains the methods for dealing with the interaction datasets from different databases and input
    """

    @staticmethod
    def biogrid_readfile(filename):
        """ Read the Biogrid tab delimited file
        Args:
            filename (:obj: str): the file path
        Returns:
            :obj:'DataFrame from pandas'        
        """

        df_raw = pd.read_csv(filename, sep='\t', skiprows=[1, 24, 1], header=25, low_memory=False)
        return df_raw

    @staticmethod
    def to_networkx_map(dataset, interactor_a, interactor_b, selfloop='n'):
        """ Create a graph from dataset by using NetworkX package
        Args:
            dataset (:obj: DataFrame from pandas): the Dataframe with header and at least two columns
            interactor_a (:obj: str): column for interactor A
            interactor_b (:obj: str): column for interactor B
            selfloop (:obj: str, optional): create a graph with or without selfloop, it should be 'y' or 'n'
        Returns:
            :obj:'Graph from networkx': the undirected graph from dataset
        Raises:
            ValueError if the column(s) is not found in the dataset
            ValueError if the selfloop is not 'y' or 'n'
            ValueError if the graph is not successfully created   
        """

        if (interactor_a not in dataset) or (interactor_b not in dataset):
            raise ValueError('the column(s) is not found in the dataset')
        elif selfloop not in ['y', 'n']:
            raise ValueError("the selfloop should be 'y' or 'n'")
        else:
            print("Now we're generating the whole biological network map...")
            if selfloop == 'y':
                nx_map = nx.from_pandas_dataframe(dataset, interactor_a, interactor_b)
            elif selfloop == 'n':
                nx_map = nx.from_pandas_dataframe(dataset, interactor_a, interactor_b)
                nx_map.remove_edges_from(nx_map.selfloop_edges())
            print("Your whole biological network map has {} nodes".format(Parameter.num_of_node(nx_map)),
                  "and {} edges.".format(Parameter.num_of_edge(nx_map)))
            return nx_map

    @staticmethod
    def input_nodelist(filename):
        """ Read the file containing only nodes
        Args:
            filename (:obj: str): the file path

        Returns:
            :obj:'list': a list of nodes
        Raise:
            ValueError if the file is not successfully loaded        
        """

        try:
            with open(filename, 'r') as file:
                datas = csv.reader(file, delimiter="\t")
                data = []
                for row in datas:
                    data.append(row)
                data = [item for sublist in data for item in sublist]
            return data
        except:
            raise ValueError('file loading error, or not found!')

    @staticmethod
    def nodes_in_map(graph_map, nodelist):
        """ Print out the nodes in the nodelist but not in the graph, and return the nodes in the nodelist
        Args:
            graph_map (:obj: Graph from networkx): the undirected graph
            nodelist (:obj: list): a list of nodes

        Returns:
            :obj:'list': a list of nodes in the map     
        """
        map_nodes = list(graph_map.nodes())
        nodes_in_graph = [item for item in nodelist if item in map_nodes]
        nodes_not_in_graph = [item for item in nodelist if item not in map_nodes]
        print(len(nodes_in_graph), 'node(s) in your node list is/are in the map.')
        print(len(nodes_not_in_graph), 'node(s) is/are not in the map.')
        if len(nodes_not_in_graph) != 0:
            print('Following are the nodes not in the map:\n', nodes_not_in_graph)

        return nodes_in_graph


class Subnetwork(object):
    """ This class contains the methods for generating subnetworks
    """

    @staticmethod
    def get_subgraph(graph_map, size):
        """ Create a subgraph from the graph by picking up nodes randomly
        Args:
            graph_map (:obj: Graph from networkx): the undirected graph
            size (:obj: int): the size of the subgraph, it should larger than 1 and smaller or equal to the size of the graph
        Returns:
            :obj:'Graph from networkx': the undirected subgraph
        Raises:
            TypeError if the size is not an integer
            ValueError if the size is not the integer between 2 and the size of the graph
        """

        if not float(size).is_integer():
            raise TypeError("the size ({}) should be a integer".format(size))
        elif (size < 2) or (size > len(list(graph_map.nodes()))):
            raise ValueError('the size should be between 2 and the size of the graph')
        else:
            nodelist = list(graph_map.nodes())
            rdnodes = random.sample(set(nodelist), size)
            nx_subg = graph_map.subgraph(rdnodes)
            return nx_subg


class Parameter(object):
    """ This class contains the methods for calculating the properties of graphs
    """

    @staticmethod
    def num_of_edge(graph):
        """ Return the number of edges of the graph
        Args:
             (:obj: Graph from networkx)
        Returns:
            :obj:'int'     
        """
        return len(list(graph.edges()))

    @staticmethod
    def num_of_node(graph):
        """ Return the number of nodes of the graph
        Args:
             graph (:obj: Graph from networkx)
        Returns:
            :obj:'int'     
        """
        return len(list(graph.nodes()))

    @staticmethod
    def max_connected(graph):
        """ Return the largest connected pattern(subgraph) of the graph
        Args:
             graph (:obj: Graph from networkx)
        Returns:
            :obj:'Graph from networkx':    
        """
        return max(nx.connected_component_subgraphs(graph), key=len)

    @staticmethod
    def nodes_divided_by_diameter(graph):
        """ Return the nodes number divided by diameter of the graph
        Args:
             graph (:obj: Graph from networkx)
        Returns:
            :obj:'float': the nodes number divided by diameter of the graph  
        """
        return Parameter.num_of_node(graph)/nx.diameter(graph)

    @staticmethod
    def edges_divided_by_diameter(graph):
        """ Return the edges number divided by diameter of the graph
        Args:
             graph (:obj: Graph from networkx)
        Returns:
            :obj:'float': the edges number divided by diameter of the graph  
        """
        return Parameter.num_of_edge(graph)/nx.diameter(graph)

    @staticmethod
    def subnets_intersecion(nodelist_1, nodelist_2):
        """ Return the intersection between two node lists
        Args:
             nodelist_1 (:obj: list): a list of nodes
             nodelist_2 (:obj: list): a list of nodes
        Returns:
            :obj:'list': the list of intersection nodes
        """
        return list(set(nodelist_1) & set(nodelist_2))

    @staticmethod
    def edge_cross(nodelist_1, nodelist_2, graph_map):
        """ Return the edges number between two subnets formed by two node lists
        Args:
             nodelist_1 (:obj: list): a list of nodes
             nodelist_2 (:obj: list): a list of nodes
             graph_map (:obj: Graph from networkx): the undirected graph
        Returns:
            :obj:'int': the edge numbers between two subnets
        """
        total_nodes = nodelist_1 + nodelist_2
        nodes_in_both = Parameter.subnets_intersecion(nodelist_1, nodelist_2)
        if not nodes_in_both:
            total_subnet = graph_map.subgraph(total_nodes)
            subnet_nodelist_1 = graph_map.subgraph(nodelist_1)
            subnet_nodelist_2 = graph_map.subgraph(nodelist_2)
            edge_between = Parameter.num_of_edge(total_subnet) - Parameter.num_of_edge(subnet_nodelist_1) -\
                           Parameter.num_of_edge(subnet_nodelist_2)
            return edge_between
        else:
            nodes_only_in_list_1 = list(set(nodelist_1) - set(nodes_in_both))
            nodes_only_in_list_2 = list(set(nodelist_2) - set(nodes_in_both))
            edge_between = Parameter.edge_cross(nodes_only_in_list_1, nodes_in_both, graph_map) +\
                           Parameter.edge_cross(nodes_only_in_list_2, nodes_in_both, graph_map) +\
                           Parameter.edge_cross(nodes_only_in_list_1, nodes_only_in_list_2, graph_map) +\
                           Parameter.num_of_edge(graph_map.subgraph(nodes_in_both))
            return edge_between

    @staticmethod
    def jaccard_node_similarity(nodelist_1, nodelist_2, graph_map):
        total_nodes = nodelist_1 + nodelist_2
        nodes_no_duplicate = list(set(total_nodes))
        nodes_in_both = list(set(nodelist_1) & set(nodelist_2))
        jaccard_node = len(nodes_in_both)/len(nodes_no_duplicate)
        return jaccard_node

    @staticmethod
    def jaccard_edge_similarity(nodelist_1, nodelist_2, graph_map):
        nodes_in_both = list(set(nodelist_1) & set(nodelist_2))
        edge_both = Parameter.num_of_edge(graph_map.subgraph(nodes_in_both))
        edge_1 = Parameter.num_of_edge(graph_map.subgraph(nodelist_1))
        edge_2 = Parameter.num_of_edge(graph_map.subgraph(nodelist_2))
        jaccard_edge = edge_both/(edge_1 + edge_2 - edge_both)
        return jaccard_edge

    @staticmethod
    def parameter_list(graphs, function):
        """ Return a list of values(or objects, based on the function) from each of multiple graphs 
        Args:
             graphs (:obj: list of Graph, Graph is an object from networkx)
             function (:any function that could return a property from the Graph)
        Returns:
            :obj:'list of values(or objects, based on the function)':    
        """
        values = list(map(lambda i: function(i), graphs))
        return values


class Stats(object):
    """ This class contains the methods for statistical analysis
    """

    @staticmethod
    def connected_pathway_sig(pathway_nodes, subnet_nodes, graph_map, num_of_sub):
        """ Give the information of how significant the supgraph connected with the pathway
        Args:
            pathway_nodes (:obj: list): a list of nodes of the pathway
            subnet_nodes (:obj: list): a list of nodes you are interested in
            graph_map (:obj: Graph from networkx): the graph of the map
            num_of_sub (:obj: int): the number of subgraphs randomly generated
        """
        print("Now we're checking the pathway nodes in the map...")
        pathway_nodes_clean = Dataclean.nodes_in_map(graph_map, pathway_nodes)
        print("Now we're checking your nodes in the map...")
        subnet_nodes_clean = Dataclean.nodes_in_map(graph_map, subnet_nodes)

        # Report the nodes are both in the pathway and your subgraph
        nodes_in_both = Parameter.subnets_intersecion(pathway_nodes_clean, subnet_nodes_clean)

        if nodes_in_both:
            print('Following are your nodes in the pathway:', nodes_in_both)

        # Report how many edges between the pathway and your subgraph
            edges_between = Parameter.edge_cross(pathway_nodes_clean, subnet_nodes_clean, graph_map)
            print('There are {} edges between the pathway and your subnet.'.format(edges_between))

            size = len(subnet_nodes_clean) - len(nodes_in_both)

        else:
            print('There are no nodes you input are in the pathway')
            edges_between = Parameter.edge_cross(pathway_nodes_clean, subnet_nodes_clean, graph_map)
            print('There are {} edges between the pathway and your subnet.'.format(edges_between))
            size = len(subnet_nodes_clean)

        nodes_in_map = list(graph_map.nodes())
        nodes_not_in_pathway = [item for item in nodes_in_map if item not in pathway_nodes_clean]
        graph_for_random = graph_map.subgraph(nodes_not_in_pathway)

        count = 0

        # Random subgraphs (exclude the nodes in the pathway)

        for i in range(num_of_sub):
            random_subnet = Subnetwork.get_subgraph(graph_for_random, size)
            random_subnet_nodelist = list(random_subnet.nodes()) + nodes_in_both
            random_edge_between = Parameter.edge_cross(pathway_nodes_clean, random_subnet_nodelist, graph_map)
            if random_edge_between >= edges_between:
                count += 1

        print('{} of'.format(count), num_of_sub, 'random subnets have more edges connected with the pathway.')
        if count == 0:
            print('Probability that chance alone gave us the number of edges connected with the pathway '
                  'at least', edges_between, 'is smaller than {0:e}'.format(1/float(num_of_sub)))
        else:
            print('Probability that chance alone gave us the number of edges connected with the pathway '
                  'at least', edges_between, 'is {0:e}'.format(count/float(num_of_sub)))

    @staticmethod
    def connected_pathway_distribution(pathway_nodes, graph_map, size, num_of_sub, intersect=[]):
        print("Now we're checking the pathway nodes in the map...")
        pathway_nodes_clean = Dataclean.nodes_in_map(graph_map, pathway_nodes)
        nodes_in_map = list(graph_map.nodes())
        nodes_not_in_pathway = [item for item in nodes_in_map if item not in pathway_nodes_clean]
        graph_for_random = graph_map.subgraph(nodes_not_in_pathway)

        edges_between = []
        for i in range(num_of_sub):
            random_subnet = Subnetwork.get_subgraph(graph_for_random, size)
            random_subnet_nodes = list(random_subnet.nodes()) + intersect
            random_subnet_edge_between = Parameter.edge_cross(random_subnet_nodes,
                                                              pathway_nodes_clean, graph_map)
            edges_between.append(random_subnet_edge_between)

        edges_between = np.array(edges_between)

        filename = 'connected_pathway_distribution_' + str(size) + '_' + str(num_of_sub) + '.txt'
        text_file = open(filename, "w")
        for i in edges_between:
            text_file.write(str(i))
            text_file.write('\n')
        text_file.close()

        weights = np.ones_like(edges_between) / float(len(edges_between))

        plt.hist(edges_between, facecolor='g', weights=weights,
                 bins=np.arange(min(edges_between), max(edges_between) + 1, 1))
        plt.title('Histogram of The Edge Number between \n The Pathway and Random Subgraphs \n'
                  '(Size: {}, '.format(size) + 'Number of Subgraphs: {})\n'.format(num_of_sub) +
                  '(Intersection: {})'.format(intersect))

        plt.xlabel('Edge Number between The Pathway and Random Subgraphs')
        plt.ylabel('Probability')
        xint = range(int(min(edges_between) - 1), int(math.ceil(max(edges_between)) + 2))
        plt.xticks(xint)
        figname = 'connected_pathway_distribution_' + str(size) + '_' + str(num_of_sub) + '.png'
        plt.savefig(figname)
        plt.close()

    @staticmethod
    def max_connected_sig(subnet_nodes, graph_map, num_of_sub):
        """ Give the information of how significant the largest connected pattern(subgraph) in your graph
        Args:
            subnet_nodes (:obj: list): a list of nodes you are interested in
            graph_map (:obj: Graph from networkx): the graph of the map
            num_of_sub (:obj: int): the number of subgraphs randomly generated
        """
        print("Now we're checking your nodes in the map...")
        subnet_nodes_clean = Dataclean.nodes_in_map(graph_map, subnet_nodes)

        subnet_nodes_clean_graph = graph_map.subgraph(subnet_nodes_clean)
        largest_of_graph = Parameter.max_connected(subnet_nodes_clean_graph)
        size_of_largest = Parameter.num_of_node(largest_of_graph)
        print('The largest connected subnet in your graph has {} nodes.'.format(size_of_largest))

        n = len(subnet_nodes_clean)
        count = 0

        print('Generating random subnets with size {}...'.format(n))
        print('Finding the largest connected pattern in each of these subnets...')
        for i in range(num_of_sub):
            random_subnet = Subnetwork.get_subgraph(graph_map, n)
            random_subnet_largest = Parameter.max_connected(random_subnet)
            random_subnet_largest_size = Parameter.num_of_node(random_subnet_largest)
            if random_subnet_largest_size >= size_of_largest:
                count += 1

        print('{} of'.format(count), num_of_sub, 'random subnets have the largest connected subgraph '
              'larger than your graph.')
        if count == 0:
            print('Probability that chance alone gave us the size of '
                  'largest connected subnet at least', size_of_largest, 'nodes is smaller than',
                  '{0:e}'.format(1/float(num_of_sub)))
        else:
            print('Probability that chance alone gave us the size of '
                  'largest connected subnet at least', size_of_largest, 'nodes is {0:e}'.format(count/float(num_of_sub)))

    @staticmethod
    def max_connected_distribution(graph_map, size, num_of_sub):
        print('Generating random subnets with size {}...'.format(size))
        print('Finding the largest connected pattern in each of these subnets...')
        largest_size = []
        for i in range(num_of_sub):
            random_subnet = Subnetwork.get_subgraph(graph_map, size)
            random_subnet_largest = Parameter.max_connected(random_subnet)
            random_subnet_largest_size = Parameter.num_of_node(random_subnet_largest)
            largest_size.append(random_subnet_largest_size)

        largest_size = np.array(largest_size)

        filename = 'max_connected_distribution_' + str(size) + '_' + str(num_of_sub) + '.txt'
        text_file = open(filename, "w")
        for i in largest_size:
            text_file.write(str(i))
            text_file.write('\n')
        text_file.close()

        weights = np.ones_like(largest_size) / float(len(largest_size))

        plt.hist(largest_size, facecolor='g', weights=weights,
                 bins=np.arange(min(largest_size), max(largest_size) + 1, 1))
        plt.title('Histogram of The Node Number of The Largest Connected Subgraphs \n'
                  '(Size: {}, '.format(size) + 'Number of Subgraphs: {})'.format(num_of_sub))
        plt.xlabel('Node Number of The Largest Connected Subgraphs')
        plt.ylabel('Probability')
        xint = range(int(min(largest_size) - 1), int(math.ceil(max(largest_size)) + 2))
        plt.xticks(xint)
        figname = 'max_connected_distribution_' + str(size) + '_' + str(num_of_sub) + '.png'
        plt.savefig(figname)
        plt.close()

    @staticmethod
    def diff_size_max_connected_sig(graph, graph_map, size, num_of_sub, num_shuffle):
        """ Give the information of how significant the largest connected pattern(subgraph) in your graph
        Args:
            graph (:obj: Graph from networkx): the graph you are interested in
            graph_map (:obj: Graph from networkx): the graph of the map
            size (:obj: int): the size of subnetworks randomly generated
            num_of_sub (:obj: int): the number of subgraphs randomly generated
            num_shuffle (:obj: int): the number of shuffle in the test
        """
        random_subnet_target = Subnetwork.get_subgraphs(graph, size, num_of_sub)
        random_subnet_map = Subnetwork.get_subgraphs(graph_map, size, num_of_sub)

        random_subnet_target_largest = Parameter.parameter_list(random_subnet_target, Parameter.max_connected)
        random_subnet_map_largest = Parameter.parameter_list(random_subnet_map, Parameter.max_connected)

        random_subnet_target_largest_size = Parameter.parameter_list(random_subnet_target_largest,
                                                                     Parameter.num_of_node)
        random_subnet_map_largest_size = Parameter.parameter_list(random_subnet_map_largest,
                                                                  Parameter.num_of_node)

        diff_mean = (sum(random_subnet_target_largest_size)/len(random_subnet_target_largest_size)) - \
                    (sum(random_subnet_map_largest_size)/len(random_subnet_map_largest_size))

        data_merge = random_subnet_target_largest_size + random_subnet_map_largest_size

        count = 0
        for i in range(0, num_shuffle):
            random.shuffle(data_merge)
            grpa = data_merge[:num_of_sub]
            grpb = data_merge[num_of_sub:]
            diff_mean_shuffle = (sum(grpa)/len(grpa)) - (sum(grpb)/len(grpb))
            if abs(diff_mean_shuffle) >= diff_mean:
                count += 1

        print("Observed difference of two means: {}".format(diff_mean))
        print(count, "out of", num_shuffle, "experiments had a difference of two means",
              "greater than or equal to", diff_mean, ".")
        print("The chance of getting a difference of two means grater than or equal to",
              diff_mean, "is", (count / float(num_shuffle)), ".")

    @staticmethod
    def percentage(values):
        """ Calculate the percentage of every values in a list
        Args:
            values (:obj: list)
        Returns:
            :obj:'list of tuples': each tuple contains (value, the percentage of this value in the list)       
        """
        c = Counter(values)
        percent = sorted([(i, c[i] / len(values) * 100.0) for i in c])
        return percent

    @staticmethod
    def get_large_graph(graph, size, number, show='n'):
        """ Get the largest graph in the graphs
        Args:
            graph (:obj: list of Graph, Graph is an object from networkx)
            size (:obj: int): the size of subgraphs, it should larger than 1 and smaller or equal to the size of the graph
            number (:obj: int): the number of subgraphs, it should it should larger than 1
            show (:obj: str): 'y' to show the graph
        Returns:
            :obj:'Graph from networkx': the largest graph in the graphs     
        """
        sub_nets = Subnetwork.get_subgraphs(graph, size, number)
        print("Now we've generated", number, 'subnetworks with', size, 'nodes.')

        graphs = Parameter.parameter_list(sub_nets, Parameter.max_connected)
        print("Now we've got the largest patterns of each of these subnetworks.")

        largest = max(graphs, key=len)
        print('The size of the largest graph in these graphs is', len(list(largest.nodes())), '.')
        if show == 'y':
            nx.draw(largest, with_labels=True)
            plt.show()
        return largest
