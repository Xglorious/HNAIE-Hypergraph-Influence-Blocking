from collections import deque
from hypergraph import Hypergraph_tool

class HNAIE:
    def __init__(self, vertices, SN, k, hyperedges, vertice_to_hedge_weight, alpha, beta, weight_dic):
        self.vertices = vertices
        self.SN = SN
        self.k = k
        self.hyperedges = hyperedges
        self.vertice_to_hedge_weight = vertice_to_hedge_weight
        self.alpha = alpha
        self.beta = beta
        self.weight_dic = weight_dic
        pass

    def CMIA_O(self):
        SP = set()  
        NegS = set() 
        DPADV = {v: 0 for v in self.vertices} 
        miia = {v: set() for v in self.vertices}
        mioa = {v: set() for v in self.vertices}

        path = {v: list() for v in self.vertices}

        for v in self.vertices:
            temp, path[v] = self.bfs_shortest_path_with_predecessors_miia(v, 1)
            miia[v] = temp.keys()
            tempath_v, paths_v = self.bfs_shortest_path_with_predecessors_mioa(v, 1)
            mioa[v] = tempath_v.keys()

        for u in self.SN:
            NegS |= mioa[u]
        NegS = NegS - self.SN

        for u in NegS:
            temp_u = self.aPN(u, SP, miia[u], path[u], 1, self.vertice_to_hedge_weight)
            for v in miia[u] - self.SN:
                DPADV[v] += temp_u - self.aPN(u, SP.union({v}), miia[u], path[u], 1, self.vertice_to_hedge_weight)

        eligible_nodes = [v for v in self.vertices if v not in self.SN.union(SP) and DPADV[v] != 0]
        SP = sorted(eligible_nodes, key=lambda v: DPADV[v], reverse=True)[:self.k]
        
        return SP


    def aPN(self, u, SP, miia_u, path, T, vertice_to_hedge_weight):
        ZP = {t: set() for t in range(T + 1)}  
        ZN = {t: set() for t in range(T + 1)}
        ZP[0] = SP.intersection(miia_u)  
        ZN[0] = self.SN.intersection(miia_u)

        pP = {v: {t: 0 for t in range(T + 1)} for v in miia_u}  
        pN = {v: {t: 0 for t in range(T + 1)} for v in miia_u}
        apP = {v: {t: 0 for t in range(T + 1)} for v in miia_u}
        apN = {v: {t: 0 for t in range(T + 1)} for v in miia_u}  
        for v in SP.intersection(miia_u):
            pP[v][0] = 1
            for t in range(T + 1):
                apP[v][t] = 1
        for v in self.SN.intersection(miia_u):
            pN[v][0] = 1
            for t in range(T + 1):
                apN[v][t] = 1
        t = 0
        while ZN[t] and t < T:
            temp_P = {v: 1 for v in miia_u}  
            temp_N = {v: 1 for v in miia_u}
            tool = Hypergraph_tool(self.hyperedges, self.vertices)
            for v in ZP[t]:
                if len(path[v]) > 1:  
                    w = path[v][1]
                    temp_P[w] *= (1 - pP[v][t] * tool.find_probability_between_nodes(v, w, self.alpha, self.vertice_to_hedge_weight))
                    ZP[t + 1].add(w)
            for v in ZP[t + 1]:
                pP[v][t + 1] = (1 - temp_P[v]) * (1 - apN[v][t]) * (1 - apP[v][t])
                for tao in range(t + 1, T + 1):
                    apP[v][tao] = apP[v][t] + pP[v][t + 1]
            for v in ZN[t]:
                if len(path[v]) > 1:  
                    w = path[v][1]
                    temp_N[w] *= (1 - pN[v][t] * tool.find_probability_between_nodes(v, w, self.beta, self.vertice_to_hedge_weight))
                    ZN[t + 1].add(w)
            for v in ZN[t + 1]:
                pN[v][t + 1] = temp_P[v] * (1 - temp_N[v]) * (1 - apN[v][t]) * (1 - apP[v][t])
                for tao in range(t + 1, T + 1):
                    apN[v][tao] = apN[v][t] + pN[v][t + 1]
                    if apN[v][tao] > 1:
                        print("WARING!!")

            t = t + 1
        return apN[u][T]   
    

    def reconstruct_path_miia(self, predecessors, start_vertex, target_vertex):
        path = [target_vertex]
        while target_vertex != start_vertex and target_vertex in predecessors:
            target_vertex = predecessors[target_vertex]
            path.append(target_vertex)
        return path if target_vertex == start_vertex else None


    def bfs_shortest_path_with_predecessors_miia(self, start_vertex, T):
        visited = {}
        queue = deque([(start_vertex, 0)])
        visited[start_vertex] = 0
        predecessors = {start_vertex: None}

        while queue:
            current_vertex, current_distance = queue.popleft()
            if current_distance > T:
                break  

            for hyperedge_key, hyperedge in self.hyperedges.items():
                if current_vertex in hyperedge:
                    for neighbor in hyperedge:
                        if neighbor == current_vertex:
                            continue 

                        if neighbor in self.SN:
                            
                            if neighbor not in visited:
                                visited[neighbor] = current_distance + 1
                                predecessors[neighbor] = current_vertex
                        else:
                            if neighbor not in visited or visited[neighbor] > current_distance + 1:
                                visited[neighbor] = current_distance + 1
                                predecessors[neighbor] = current_vertex
                                queue.append((neighbor, current_distance + 1))

        reachable_vertices = {vertex: dist for vertex, dist in visited.items() if dist <= T}
        paths = {vertex: self.reconstruct_path_miia(predecessors, start_vertex, vertex)
                for vertex in reachable_vertices if vertex in predecessors}

        return reachable_vertices, paths

    def reconstruct_path_mioa(self, predecessors, start_vertex, target_vertex):
        path = [target_vertex]  
        while target_vertex != start_vertex and target_vertex in predecessors:
            target_vertex = predecessors[target_vertex]  
            if target_vertex is not None:
                path.append(target_vertex)  
        if target_vertex == start_vertex:
            return path
        else:
            return None  


    def bfs_shortest_path_with_predecessors_mioa(self, start_vertex, T):
        visited = {}  
        queue = deque([(start_vertex, 0)])  
        visited[start_vertex] = 0  
        predecessors = {start_vertex: None}  

        while queue:
            current_vertex, current_distance = queue.popleft()
            if current_distance > T:
                continue  

            for hyperedge_key, hyperedge in self.hyperedges.items():
                if current_vertex in hyperedge:
                    for neighbor in hyperedge:
                        if neighbor != current_vertex and (
                                neighbor not in visited or visited[neighbor] > current_distance + 1):
                            visited[neighbor] = current_distance + 1
                            predecessors[neighbor] = current_vertex
                            queue.append((neighbor, current_distance + 1))

            reachable_vertices = {vertex: visited[vertex] for vertex, distance in visited.items() if distance <= T}

            paths = {}
            for vertex in reachable_vertices:
                path = self.reconstruct_path_mioa(predecessors, start_vertex, vertex)
                if path:  
                    paths[vertex] = path

        return reachable_vertices, paths