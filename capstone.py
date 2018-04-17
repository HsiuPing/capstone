# Project: BMI Capstone Project
# Filename: capstone.py
# Time: 201801
# Editor: Hsiu-Ping Lin

import csv
import random
import numpy as np
import math
import mygene
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
            if selfloop == 'y':
                nx_map = nx.from_pandas_dataframe(dataset, interactor_a, interactor_b)
            elif selfloop == 'n':
                nx_map = nx.from_pandas_dataframe(dataset, interactor_a, interactor_b)
                nx_map.remove_edges_from(nx_map.selfloop_edges())
            print("Retrieve the protein-protein interaction network with {} nodes".format(Parameter.num_of_node(nx_map)),
                  "and {} edges.".format(Parameter.num_of_edge(nx_map)))
            return nx_map

    @staticmethod
    def id_mapping(filename):
        """ Map the gene id to the gene name from the pathway file from wikipathway
        Args:
            filename (:obj: str): the file path
        Returns:
            :obj:'list': a list of gene names        
        """
        df_id = pd.read_csv(filename, sep='\t', header=0, low_memory=False)
        mg = mygene.MyGeneInfo()
        gene_list = []
        count = 0
        gene_not_found = 0
        print("Now we're mapping the gene list of the pathway.")
        for index, row in df_id.iterrows():
            gene_id = row['Identifier']
            gene_name = mg.getgene(gene_id, fields='symbol')
            if gene_name is not None:
                gene_list.append(gene_name['symbol'])
            else:
                print('ID: ', gene_id, 'is not found.')
                gene_not_found += 1
            count += 1

        print("There is/are", count, "gene(s) in this pathway.", gene_not_found,
              "of them is/are not found when id mapping.")
        return gene_list

    @staticmethod
    def input_nodelist(filename, printout=False):
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
                if printout is True:
                    print("Loading the file:", filename)
                data = []
                for row in datas:
                    data.append(row)
                data = [item for sublist in data for item in sublist]
            return data
        except:
            raise ValueError('file loading error, or not found!')

    @staticmethod
    def input_distribution(filename):
        """ Read the file containing the distribution values
                Args:
                    filename (:obj: str): the file path

                Returns:
                    :obj:'list': a list of float representing the distribution
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
                data = [float(i) for i in data]
            return data
        except:
            raise ValueError('file loading error, or not found!')

    @staticmethod
    def nodes_in_map(graph_map, nodelist, printout=False):
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
        if printout is True:
            print(len(nodes_in_graph), 'genes in your gene set are in the PPI network.')
            print(len(nodes_not_in_graph), 'genes are not in the PPI network.')
        if len(nodes_not_in_graph) != 0:
            if printout is True:
                print('Following are the genes not in the PPI network:', nodes_not_in_graph)

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
    def edge_cross_normalized(nodelist_1, nodelist_2, graph_map):
        """ Return the value normalized by number of nodes in nodelist_1
        Args:
             nodelist_1 (:obj: list): a list of nodes
             nodelist_2 (:obj: list): a list of nodes
             graph_map (:obj: Graph from networkx): the undirected graph
        Returns:
            :obj:'int': the edge numbers between two subnets
        """
        edge_between = Parameter.edge_cross(nodelist_1, nodelist_2, graph_map)
        edge_between_normalized = edge_between/float(len(nodelist_1))
        return edge_between_normalized

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
        """ Give the information of how significant the subgraph connected with the pathway
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
    def connected_pathway_distribution(pathway_nodes, graph_map, num_of_sub, size=None, pathway_name=None):
        """ Generate the text file and png file to show the distribution of the edge number between the given pathway
            and the random subgraphs
        Args:
            pathway_nodes (:obj: list): a list of nodes of the pathway
            graph_map (:obj: Graph from networkx): the graph of the map
            num_of_sub (:obj: int): the number of subgraphs randomly generated
            size (:obj: int or None, optional): the size of the subgraph, it should larger than 1 and smaller or equal
                                                to the size of the graph. If it is None, subgraphs with size from 2 to 999
                                                will be randomly generated.
            pathway_name (:obj: str or None, optional): the name of your pathway
            
        """
        print("Now we're checking the pathway nodes in the map...")
        pathway_nodes_clean = Dataclean.nodes_in_map(graph_map, pathway_nodes)

        edges_between = []

        if size is None:
            for i in range(num_of_sub):
                size = random.choice(list(range(2, 1000)))
                random_subnet = Subnetwork.get_subgraph(graph_map, size)
                random_subnet_nodes = list(random_subnet.nodes())
                random_subnet_edge_between = Parameter.edge_cross_normalized(random_subnet_nodes,
                                                                             pathway_nodes_clean, graph_map)

                edges_between.append(random_subnet_edge_between)

            edges_between = np.array(edges_between)

            filename = 'connected_pathway_distribution_{}_random_size_'.format(pathway_name) + str(num_of_sub) + '.txt'

        else:
            for i in range(num_of_sub):
                random_subnet = Subnetwork.get_subgraph(graph_map, size)
                random_subnet_nodes = list(random_subnet.nodes())
                random_subnet_edge_between = Parameter.edge_cross_normalized(random_subnet_nodes,
                                                                             pathway_nodes_clean, graph_map)

                edges_between.append(random_subnet_edge_between)

            edges_between = np.array(edges_between)

            filename = 'connected_pathway_distribution_{}_{}_'.format(pathway_name, size) + str(num_of_sub) + '.txt'

        text_file = open(filename, "w")
        for i in edges_between:
            text_file.write(str(i))
            text_file.write('\n')
        text_file.close()

        return edges_between

    @staticmethod
    def connected_pathway_distribution_figure(distribution, size='random', name=None):
        # Get current size
        fig_size = plt.rcParams["figure.figsize"]

        # Set figure width to 12 and height to 9
        fig_size[0] = 12
        fig_size[1] = 9
        plt.rcParams["figure.figsize"] = fig_size
        plt.rcParams["patch.force_edgecolor"] = True

        weights = np.ones_like(distribution) / float(len(distribution))

        plt.hist(distribution, facecolor='g', weights=weights,
                 bins=np.arange(min(distribution), math.ceil(max(distribution)) + 0.1, 0.025))
        title = 'Histogram of The Normalized Edge Number between ' + str(name) + ' and Random Subgraphs \n' \
                + 'Number of Subgraphs:' + str(len(distribution))
        plt.title(title)

        if math.ceil(max(distribution)) > 10:
            xint = np.arange(min(distribution), math.ceil(max(distribution)) + 0.1, 0.5)
        elif (math.ceil(max(distribution)) <= 10) and (math.ceil(max(distribution)) > 5):
            xint = np.arange(min(distribution), math.ceil(max(distribution)) + 0.1, 0.2)
        else:
            xint = np.arange(min(distribution), math.ceil(max(distribution)) + 0.1, 0.1)
        plt.xticks(xint)

        plt.xlabel('Edge Number between The Pathway and Random Subgraphs\n Normalized by size of each subgraph')
        plt.ylabel('Probability')
        plt.xticks(rotation=90)
        figname = 'connected_pathway_distribution_' + str(name) + '_' + str(size) + '_' +\
                  str(len(distribution)) + '.png'
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
        """ Generate the text file to show the distribution of the node number of the largest linked
            subgraph from each of the random subgraphs
        Args:
            graph_map (:obj: Graph from networkx): the graph of the map
            size (:obj: int): the size of the subgraph, it should larger than 1 and smaller or equal to the size of the graph
            num_of_sub (:obj: int): the number of subgraphs randomly generated
        """
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
        return largest_size

    @staticmethod
    def max_connected_distribution_figure(distribution, size, value=None):
        # Get current size
        fig_size = plt.rcParams["figure.figsize"]

        # Set figure width to 12 and height to 9
        fig_size[0] = 18
        fig_size[1] = 9
        plt.rcParams["figure.figsize"] = fig_size
        plt.rcParams["patch.force_edgecolor"] = True
        if value is not None:
            distribution_with_value = distribution + [value]
        else:
            distribution_with_value = distribution

        weights = np.ones_like(distribution) / float(len(distribution))

        plt.hist(distribution, facecolor='g', weights=weights,
                 bins=np.arange(min(distribution), max(distribution) + 1, 1))
        if value is not None:
            plt.axvline(value, color='r', linestyle='dashed', linewidth=3)

        plt.title('Histogram of The Node Number of The Largest Connected Subgraphs \n'
                  '(Size: {}, '.format(size) + 'Number of Subgraphs: {})'.format(len(distribution)), fontsize=28)
        plt.xlabel('Node Number of The Largest Connected Subgraphs', fontsize=26)
        plt.ylabel('Probability', fontsize=26)
        xint = range(int(min(distribution_with_value) - 1), int(math.ceil(max(distribution_with_value)) + 2))
        plt.xticks(xint, fontsize=12)
        plt.yticks(fontsize=16)
        plt.xticks(rotation=90)
        figname = 'max_connected_distribution_' + str(size) + '_' + str(len(distribution)) + '.png'
        plt.savefig(figname)
        plt.close()

    @staticmethod
    def p_in_distribution(value, distribution):
        count = 0
        for i in distribution:
            if i >= value:
                count += 1

        print('{} of'.format(count), len(distribution), 'random subnets have the value '
              'larger than your graph.')
        if count == 0:
            print('Probability that chance alone gave us that larger or equal to your value', value,
                  'in this distribution is smaller than',
                  '{0:e}'.format(1/float(len(distribution))))
            return 0
        else:
            print('Probability that chance alone gave us that larger or equal to your value', value,
                  'in this distribution is {0:e}'.format(count/float(len(distribution))))
            return count/float(len(distribution))

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
