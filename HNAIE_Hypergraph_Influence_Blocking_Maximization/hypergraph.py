
class Hypergraph_tool:
    def __init__(self, hyperedges, vertices):
        self.hyperedges = hyperedges
        self.vertices = vertices
        pass

    def compute_degrees(self):
        degrees = {}
        for hyperedge in self.hyperedges.values():
            for node in hyperedge:
                if node not in degrees:
                    degrees[node] = set()
                degrees[node].update(hyperedge)
        for node in degrees:
            degrees[node] = len(degrees[node]) - 1
        return degrees
    
    def get_top_k_nodes_by_degree(self, k):
        degrees = self.compute_degrees()
        sorted_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
        top_k_nodes = [node[0] for node in sorted_nodes[:k]]
        return top_k_nodes
    
    def compute_hyperdegrees(self):
        hyperdegrees = {}
        for hyperedge in self.hyperedges.values():
            for node in hyperedge:
                if node in hyperdegrees:
                    hyperdegrees[node] += 1
                else:
                    hyperdegrees[node] = 1
        return hyperdegrees

    def get_top_k_nodes(self, k):
        hyperdegrees = self.compute_hyperdegrees()
        sorted_nodes = sorted(hyperdegrees.items(), key=lambda x: x[1], reverse=True)
        top_k_nodes = [node[0] for node in sorted_nodes[:k]]
        return top_k_nodes
    
    def find_hyperedges_containing_node(self, u):
        hyperedges_containing_u = {key for key, nodes in self.hyperedges.items() if u in nodes}
        return hyperedges_containing_u

    def find_neighbors(self, u):
        u_neighbors = set()
        u_hyperedge = self.find_hyperedges_containing_node(u)
        for Hedge in u_hyperedge:
            u_neighbors |= set(self.hyperedges[Hedge])
        return u_neighbors
    
    def find_probability_between_nodes(self, u, v, theta, vertice_to_hedge_weight):
        common_hyperedges = [he for he in self.hyperedges if u in self.hyperedges[he] and v in self.hyperedges[he]]
        probability = 1
        for edge in common_hyperedges:
            probability_u_to_edge = vertice_to_hedge_weight[u][edge]
            probability *= (1 - probability_u_to_edge * theta)

        return 1 - probability

    def cal_temp_weight_key(self, u, weight_dic):
        u_hyperedges = self.find_hyperedges_containing_node(u)
        weight_key = {hyedge: weight_dic[hyedge] for hyedge in u_hyperedges}
        return weight_key
    
    