import random
from hypergraph import Hypergraph_tool

class HCIC_model:
    def __init__(self, SN, SP, vertices, hyperedges, alpha, beta):
        self.SN = SN
        self.SP = SP
        self.vertices = vertices
        self.hyperedges = hyperedges
        self.alpha = alpha
        self.beta = beta
    
    def extend_hyper_IC(self, vertice_to_hedge_weight):
        Negative_nodes = set(self.SN)
        Positive_nodes = set(self.SP)
        Inactive_nodes = set(self.vertices) - Negative_nodes - Positive_nodes
        Negative_nodes_number = [len(Negative_nodes)]
        Positive_nodes_number = [len(Positive_nodes)]
        this_time_Negative_nodes = Negative_nodes
        this_time_Positive_nodes = Positive_nodes

        while this_time_Positive_nodes or this_time_Negative_nodes:
            new_Negative_nodes = set()
            new_Positive_nodes = set()
            tool = Hypergraph_tool(self.hyperedges, self.vertices)

            for p_node in list(this_time_Positive_nodes):
                for i_edge in tool.find_hyperedges_containing_node(p_node):
                    if random.random() < vertice_to_hedge_weight[p_node][i_edge]:  
                        for i_node in self.hyperedges[i_edge].intersection(Inactive_nodes):
                            if random.random() < self.alpha:  
                                new_Positive_nodes.add(i_node)
                                Positive_nodes.add(i_node)
                                Inactive_nodes -= {i_node}
            for n_node in list(this_time_Negative_nodes):
                for i_edge in tool.find_hyperedges_containing_node(n_node):
                    if random.random() < vertice_to_hedge_weight[n_node][i_edge]:  
                        for i_node in self.hyperedges[i_edge].intersection(Inactive_nodes):
                            if random.random() < self.beta:  
                                new_Negative_nodes.add(i_node)
                                Negative_nodes.add(i_node)
                                Inactive_nodes -= {i_node}

            Positive_nodes_number.append(len(new_Positive_nodes))
            Negative_nodes_number.append(len(new_Negative_nodes))
            this_time_Positive_nodes = new_Positive_nodes
            this_time_Negative_nodes = new_Negative_nodes

        return Positive_nodes_number, Negative_nodes_number
