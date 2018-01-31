# Project: BMI Capstone Project
# Filename: capstone.py
# Time: 201801
# Editor: Hsiu-Ping Lin

import csv
import random
from scipy import stats
from matplotlib import pyplot as plt
import pandas as pd
from collections import Counter
import networkx as nx


# APG: Standards compliance, Robustness: Python 3 standard declares classes with
#   class ClassName(object), to make it explicit that they derive from object
class Dataclean:
    """ This class contains the methods for dealing with the interaction datasets from different databases and input
    """

    @staticmethod
    def biogrid_readfile(filename):
        """ Read the Biogrid tab delimited file
        Args:
            filename (:obj: str): the file path

        Returns:
            :obj:'DataFrame from pandas'
        Raise:
            ValueError if the file is not successfully loaded
        """
        try:
            df_raw = pd.read_csv(filename, sep='\t', skiprows=[1, 24, 1], header=25, low_memory=False)
            return df_raw
        except:
            # APG: Robustness: include more information in the error message, such as filename, type of error
            # in fact, since pandas has robust error reporting, it would be easiest to not catch the
            # exception, and simply let the pandas error be raised in the calling code
            # alternatively, you could reraise the pandas error
            # try
            #    pd.read_csv()
            # except Exception as e:
            #    raise ValueError("unable to read file" + e)
            raise ValueError('file loading error, or not found!')

    @staticmethod
    def tonetworkx(dataset, colname_a, colname_b, selfloop='n'):
        """ Create a graph from dataset by using NetworkX package
        Args:
            dataset (:obj: DataFrame from pandas): the Dataframe with header and at least two columns
            colname_a (:obj: str): column for interactor A
            colname_b (:obj: str): column for interactor B
            selfloop (:obj: str): create a graph with or without selfloop, it should be 'y' or 'n'
        Returns:
            :obj:'Graph from networkx': the undirected graph from dataset
        Raises:
            ValueError if the column(s) is not found in the dataset
            ValueError if the selfloop is not 'y' or 'n'
            ValueError if the graph is not successfully created
        """
        # APG: Documentation: to indicate an optional argument use 'optional'; e.g.:
        # selfloop (:obj: str, optional):

        if (colname_a not in dataset) or (colname_b not in dataset):
            raise ValueError('the column(s) is not found in the dataset')
        elif selfloop not in ['y', 'n']:
            raise ValueError("the selfloop should be 'y' or 'n'")
        else:
            if selfloop == 'y':
                nx_map = nx.from_pandas_edgelist(dataset, colname_a, colname_b)
            elif selfloop == 'n':
                nx_map = nx.from_pandas_edgelist(dataset, colname_a, colname_b)
                nx_map.remove_edges_from(nx_map.selfloop_edges())
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
        # APG: Documentation: small point: docstrings typically put a blank line between each section
        # e.g., see 'def fetch_bigtable_rows' example in https://google.github.io/styleguide/pyguide.html

        try:
            with open(filename, 'r') as file:
                datas = csv.reader(file, delimiter="\t")
                data = []
                for row in datas:
                    data.append(row)
                data = [item for sublist in data for item in sublist]
            return data
        except:
            # APG: Robustness: basically same Robustness as above
            raise ValueError('file loading error, or not found!')

    @staticmethod
    def nodeinmap(graph, nodelist):
        """ Print out the nodes in the nodelist but not in the graph, and return the subgraph by nodes in the nodelist
        Args:
            graph (:obj: Graph from networkx): the undirected graph
            nodelist (:obj: list): a list of nodes

        Returns:
            :obj:'Graph from networkx': the subgraph by nodes in nodelist
        Raise:
            ValueError if the file is not successfully loaded
        """
        graphnode = list(graph.nodes())
        nodesingraph = list(set(graphnode) & set(nodelist))
        nodesnotingraph = [item for item in nodelist if item not in nodesingraph]
        # APG: Readability: names like 'nodesnotingraph' are hard to read; use 'nodes_not_in_graph' instead
        # these are called 'lower_with_under' names; see https://google.github.io/styleguide/pyguide.html again
        # you can improve many names with this concept
        print(len(nodesingraph), 'nodes in your node list are in the graph.')
        print(len(nodesnotingraph), 'nodes are not in the graph.')
        print('Following are the nodes not in the graph:\n', nodesnotingraph)

        nx_subg = graph.subgraph(nodesingraph)
        return nx_subg

class Subgraph:
    """ This class contains the methods for generating subnetworks
    """

    @staticmethod
    def subgraph(graph, size):
        """ Create a subgraph from the graph by picking up nodes randomly
        Args:
            graph (:obj: Graph from networkx): the undirected graph
            size (:obj: int): the size of the subgraph, it should larger than 1 and smaller or equal to the size of the graph
        Returns:
            :obj:'Graph from networkx': the undirected subgraph
        Raises:
            ValueError if the size is not the integer between 2 and the size of the graph
            ValueError if the subgraph is not successfully created
        """

        if not float(size).is_integer():
            raise ValueError('the size should be a integer')
            # AGP: Robustness: could report the bad size with
            #   ValueError("the size ({}) should be a integer".format(size))
        elif (size < 2) or (size > len(list(graph.nodes()))):
            raise ValueError('the size should be between 2 and the size of the graph')
        else:
            try:
                nodelist = list(graph.nodes())
                rdnodes = random.sample(set(nodelist), size)
                nx_subg = graph.subgraph(rdnodes)
                return nx_subg
            except:
                # AGP: Robustness: needs more explanation, so user can solve the problem
                raise ValueError('the subgraph is not successfully created')

    @staticmethod
    def subgraph_link(graph, size):
        """ Create a linked subgraph from the graph by picking up the neighboring nodes randomly
        Args:
            graph (:obj: Graph from networkx): the undirected graph
            size (:obj: int): the size of the subgraph, it should larger than 1 and smaller or equal to the size of the graph
        Returns:
            :obj:'Graph from networkx': the undirected linked subgraph
        Raises:
            ValueError if the size is not the integer between 2 and the size of the graph
            ValueError if the subgraph is not successfully created
        """

        # APG: Functionality: these tests should check whether size > len(graph.largest_connected_subgraph())
        if not float(size).is_integer():
            raise ValueError('the size should be a integer')
        elif (size < 2) or (size > len(list(graph.nodes()))):
            raise ValueError('the size should be between 2 and the size of the graph')

        else:
            # APG: Functionality: it's not clear whether this algorithm is creating an UNBIASED
            # random connected subgraph; first, we need to define the algorithm's goal
            # it could be "choose a subgraph so all possible connected subgraphs are equally likely to be chosen"
            try:
                dic_adjlist = nx.to_dict_of_lists(graph)
                temp_k = random.choice(list(dic_adjlist))
                temp_list = []
                rdnodes = [temp_k]

                i = 1

                while i < size:
                    for node in dic_adjlist[temp_k]:
                        temp_list.append(node)

                    temp_list = [x for x in temp_list if x != temp_k]
                    temp_list = [x for x in temp_list if x not in rdnodes]

                    if temp_list:
                        temp_v = random.choice(temp_list)
                        rdnodes.append(temp_v)
                        temp_k = temp_v
                        i += 1
                    else:
                        # APG: Functionality: since the while condition is i < size but i is being reset to 1
                        # this algorithm might have a non-zero probability of being an infinite loop
                        i = 1
                        temp_k = random.choice(list(dic_adjlist))
                        temp_list = []
                        rdnodes = [temp_k]
                        continue
                nx_subg = graph.subgraph(rdnodes)
                return nx_subg
            except:
                # APG: Robustness: same comment as above, plus does this method randomly fail?
                raise ValueError('the subgraph is not successfully created')

    @staticmethod
    def multi_subgraphs(graph, size, number, linked='unlinked'):
        """ Create multiple subgraphs from the graph
        Args:
            graph (:obj: Graph from networkx): the undirected graph
            size (:obj: int): the size of subgraphs, it should larger than 1 and smaller or equal to the size of the graph
            number (:obj: int): the number of subgraphs, it should it should larger than 1
            linked (:obj: str): create linked or unlinked subgraphs, it should be 'linked' or 'unlinked', the default is unlinked
        Returns:
            :obj:'list': list of undirected subgraphs (Graph from networkx)
        Raises:
            ValueError if the size is not the integer between 2 and the size of the graph
            ValueError if the number is not the integer larger than 1
            ValueError if the subgraphs are not successfully created
        """

        if not float(size).is_integer():
            raise ValueError('the size should be a integer')
        elif not float(number).is_integer():
            raise ValueError('the number should be a integer')
        elif (size < 2) or (size > len(list(graph.nodes()))):
            raise ValueError('the size should be between 2 and the size of the graph')
        elif number < 2:
            raise ValueError('the number should be larger than 1')
        elif linked not in ['linked', 'unlinked']:
            raise ValueError('should be linked or unlinked')
        else:
            try:
                if linked == 'linked':
                    subgraphs = []
                    for i in range(number):
                        subgraphs.append(Subgraph.subgraph_link(graph, size))
                    return subgraphs
                elif linked == 'unlinked':
                    subgraphs = []
                    for i in range(number):
                        subgraphs.append(Subgraph.subgraph(graph, size))
                    return subgraphs
            except:
                raise ValueError('the subgraphs are not successfully created')


class Parameter:
    """ This class contains the methods for calculating the properties of graphs
    """

    @staticmethod
    def numofedge(graph):
        """ Return the number of edges of the graph
        Args:
             (:obj: Graph from networkx)
        Returns:
            :obj:'int'
        """
        # APG: Readability: simplify by using just one line:
        # return len(list(graph.edges()))
        value = len(list(graph.edges()))
        return value

    @staticmethod
    def numofnode(graph):
        """ Return the number of nodes of the graph
        Args:
             graph (:obj: Graph from networkx)
        Returns:
            :obj:'int'
        """
        value = len(list(graph.nodes()))
        return value

    @staticmethod
    def maxconnected(graph):
        """ Return the largest connected pattern(subgraph) of the graph
        Args:
             graph (:obj: Graph from networkx)
        Returns:
            :obj:'Graph from networkx':
        """
        maxsub = max(nx.connected_component_subgraphs(graph), key=len)
        return maxsub

    @staticmethod
    def parameterlist(graphs, function):
        """ Return a list of values(or objects, based on the function) from each of multiple graphs
        Args:
             graphs (:obj: list of Graph, Graph is an object from networkx)
             function (:any function that could return a property from the Graph)
        Returns:
            :obj:'list of values(or objects, based on the function)':
        """
        values = list(map(lambda i: function(i), graphs))
        return values


class Stats:
    """ This class contains the methods for statistical analysis
    """

    @staticmethod
    def p_valuefromz(population, data):
        """ Calculate the p-value of data based on the population and the z-score
        Args:
            population (:obj: list)
            data (:obj: float)
        Returns:
            :obj:'float': the p-value of the data based on the population
        """
        popplusdata = population + [data]
        z = stats.zscore(popplusdata)
        p_values = stats.norm.sf(abs(z))
        return p_values[-1]

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
    def getlargegraph(graphs, show='n'):
        """ Get the largest graph in the graphs
        Args:
            graphs (:obj: list of Graph, Graph is an object from networkx)
            show (:obj: str): 'y' to show the graph
        Returns:
            :obj:'Graph from networkx': the largest graph in the graphs
        """
        largest = max(graphs, key=len)
        print('The size of the largest graph is', len(list(largest.nodes())))
        if show == 'y':
            nx.draw(largest, with_labels=True)
            plt.show()
        return largest

# AGP: Readability: I note that all your methods are static, which means that these classes
# store no state; typically one would want to combine state with the methods that operate on it
# in abstract data types (ADTs)
