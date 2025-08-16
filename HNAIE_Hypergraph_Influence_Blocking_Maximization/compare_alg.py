import pandas as pd
import random
import numpy as np
import copy
import networkx as nx
import time
from hypergraph import Hypergraph_tool
from HCIC_diffusion_model import HCIC_model
import heapq
from tqdm import tqdm


class Algs:
    def __init__(self, hyperedges, vertices, dataLoader):
        self.hyperedges = hyperedges
        self.vertices = vertices
        self.dataLoader = dataLoader

        pass

    def get_top_k_nodes_excluding_sn(self, k, SN):
        tool = Hypergraph_tool(self.hyperedges, self.vertices)
        degrees = tool.compute_degrees()

        filtered_degrees = {node: degree for node, degree in degrees.items() if node not in SN}

        sorted_nodes = sorted(filtered_degrees.items(), key=lambda x: x[1], reverse=True)
        top_k_nodes = [node[0] for node in sorted_nodes[:k]]

        return top_k_nodes


    def get_top_k_hyperdegree_nodes_excluding_sn(self, k, SN):
        tool = Hypergraph_tool(self.hyperedges, self.vertices)

        hyperdegrees = tool.compute_hyperdegrees()


        filtered_hyperdegrees = {node: degree for node, degree in hyperdegrees.items() if node not in SN}


        sorted_nodes = sorted(filtered_hyperdegrees.items(), key=lambda x: x[1], reverse=True)
        top_k_nodes = [node[0] for node in sorted_nodes[:k]]

        return top_k_nodes

    def random_select_nodes(self, k, SN):
        
        all_nodes = list(self.vertices - SN)
        
        selected_nodes = random.sample(all_nodes, k)

        return selected_nodes
    

    def generateMap(self, path):
        node_dict = {}
        reverse_node_dict = {}
        df = pd.read_csv(path, index_col=False, header = None)

        arr = df.values
        node_list = []
        for each in arr:
            node_list.extend(list(map(int, each[0].split(' '))))
        node_arr = np.unique(np.array(node_list))
        for i, node in enumerate(node_arr):
            node_dict[node] = i
            reverse_node_dict[i] = node
        return node_dict, reverse_node_dict, len(node_arr), len(arr)

    def changeEdgeToMatrix(self, path):

        node_dict , reverse_node_dict,N, M = self.generateMap(path)

        matrix = np.random.randint(0, 1, size=(N, M))

        df = pd.read_csv(path, index_col=False, header=None)
        arr = df.values
        index = 0
        for each in arr:

            edge_list = list(map(int, each[0].split(" ")))
            for edge in edge_list:

                matrix[node_dict[edge]][index] = 1
            index = index + 1

        return pd.DataFrame(matrix)

    def getDegreeList(self, degree):

        matrix = []
        matrix.append(np.arange(len(degree)))
        matrix.append(degree)
        df_matrix = pd.DataFrame(matrix)
        df_matrix.index = ['node_index', 'node_degree']
        return df_matrix.sort_values(by=df_matrix.index.tolist()[1], ascending=False, axis=1)


    def getMaxDegreeNode(self, degree, seeds):

        degree_copy = copy.deepcopy(degree)
        global chosedNode
        while 1:
            flag = 0
            degree_matrix = self.getDegreeList(degree_copy)
            node_index = degree_matrix.loc['node_index']
            for node in node_index:
                if node not in seeds:
                    chosedNode = node
                    flag = 1
                    break
            if flag == 1:
                break
        return [chosedNode]

    def getTotalAdj(self, df_hyper_matrix, N):
        deg_list = []
        nodes_arr = np.arange(N)
        for node in nodes_arr:
            node_list = []
            edge_set = np.where(df_hyper_matrix.loc[node] == 1)[0]
            for edge in edge_set:
                node_list.extend(list(np.where(df_hyper_matrix[edge] == 1)[0]))
            node_set = np.unique(np.array(node_list))
            deg_list.append(len(list(node_set)) - 1)
        return np.array(deg_list)



    def updateDeg_hur(self, degree, chosenNode, df_hyper_matrix, seeds):

        edge_set = np.where(df_hyper_matrix.loc[chosenNode] == 1)[0]
        adj_set = []
        for edge in edge_set:
            adj_set.extend(list(np.where(df_hyper_matrix[edge] == 1)[0]))
        adj_set_unique = np.unique(np.array(adj_set))

        for adj in adj_set_unique:
            adj_edge_set = np.where(df_hyper_matrix.loc[adj] == 1)[0]
            adj_adj_set = []
            for each in adj_edge_set:
                adj_adj_set.extend(list(np.where(df_hyper_matrix[each] == 1)[0]))
            if adj in adj_adj_set:
                adj_adj_set.remove(adj)
            sum = 0
            for adj_adj in adj_adj_set:
                if adj_adj in seeds:
                    sum += 1

            degree[adj] -= sum

    def getSeeds_hdd(self, df_hyper_matrix,N, K, reverse_node_dict, SN):
        seeds = []
        degree = self.getTotalAdj(df_hyper_matrix, N)
        for j in range(1, K + 1):
            chosenNode = self.getMaxDegreeNode(degree, seeds)[0]

            while str(reverse_node_dict[chosenNode]) in SN:
                degree[chosenNode] = -100
                chosenNode = self.getMaxDegreeNode(degree, seeds)[0]

            seeds.append(chosenNode)
            self.updateDeg_hur(degree, chosenNode, df_hyper_matrix, seeds)

        return [str(reverse_node_dict[node]) for node in seeds]  


    #HRIS
    def getSeeds_ris(self, N, K, lamda, theta, df_hyper_matrix, reverse_node_dict, SN):
        S = []
        U = []
        SN_new = set()
        for u in SN:
            SN_new.add(int(u))
     
        for theta_iter in range(theta):
            df_matrix = copy.deepcopy(df_hyper_matrix)
 
            selected_node = random.sample(list(np.arange(len(df_hyper_matrix.index.values))), 1)[0]
   
            all_edges = np.arange(len(df_hyper_matrix.columns.values))
            prob = np.random.random(len(all_edges))
            index = np.where(prob > lamda)[0]
            for edge in index:
                df_matrix[edge] = 0

            adj_matrix = np.dot(df_matrix, df_matrix.T)
            adj_matrix[np.eye(N, dtype=np.bool_)] = 0
            df_adj_matrix = pd.DataFrame(adj_matrix)
            df_adj_matrix[df_adj_matrix > 0] = 1
            G = nx.from_numpy_array(df_adj_matrix.values)
            shortest_path = nx.shortest_path(G, source=selected_node)
            RR = []
            for each in shortest_path:
                RR.append(each)
            U.append(list(np.unique(np.array(RR))))

        for k in range(K):
            U_list = []
            for each in U:
                U_list.extend(each)
            dict_count = {}
            for each in U_list:
                if each in dict_count:
                    dict_count[each] += 1
                else:
                    dict_count[each] = 1
            
            candidate_list = sorted(dict_count.items(), key=lambda item: item[1], reverse=True)
    
            for node, _ in candidate_list:
                if str(reverse_node_dict[node]) not in SN and str(reverse_node_dict[node]) not in S:
                    S.append(str(reverse_node_dict[node]))
                    
                    break

    
            U = [each for each in U if node not in each]
        return S

    #HCI
    def computeCI(self, l, N, df_hyper_matrix):
        CI_list = []
        degree = df_hyper_matrix.sum(axis=1)
        M = len(df_hyper_matrix.columns.values)
        for i in range(0, N):
            
            edge_set = np.where(df_hyper_matrix.loc[i] == 1)[0]
            
            if l == 1:
                node_list = []
                for edge in edge_set:
                    node_list.extend(list(np.where(df_hyper_matrix[edge] == 1)[0]))
                   
                if i in node_list:
                    node_list.remove(i)
                node_set = np.unique(np.array(node_list))
                
            elif l == 2:
                node_list = []
                for edge in edge_set:
                    node_list.extend(list(np.where(df_hyper_matrix[edge] == 1)[0]))
                if i in node_list:
                    node_list.remove(i)
                node_set1 = np.unique(np.array(node_list))
                node_list2 = []
                edge_matrix = np.dot(df_hyper_matrix.T, df_hyper_matrix)
                edge_matrix[np.eye(M, dtype=np.bool_)] = 0
                df_edge_matrix = pd.DataFrame(edge_matrix)
                adj_edge_list = []
                for edge in edge_set:
                    adj_edge_list.extend(list(np.where(df_edge_matrix[edge] != 0)[0]))
                adj_edge_set = np.unique(np.array(adj_edge_list))
                for each in adj_edge_set:
                    node_list2.extend(list(np.where(df_hyper_matrix[each] == 1)[0]))
                node_set2 = list(np.unique(np.array(node_list2)))
                for node in node_set2:
                    if node in list(node_set1):
                        
                        node_set2.remove(node)
                node_set = np.array(node_set2)
            ki = degree[i]
            sum = 0
            for u in node_set:
                sum = sum + (degree[u] - 1)
            CI_i = (ki - 1) * sum
            CI_list.append(CI_i)
        return CI_list


    def getSeeds_ci(self, l, N, K, df_hyper_matrix, reverse_node_dict, SN, alpha, beta, vertice_to_hedge_weight):
        a = time.time()
        
        seeds = []
        n = np.ones(N)  
        CI_list = self.computeCI(l, N, df_hyper_matrix)
        CI_arr = np.array(CI_list)

   
        SN_set = set(SN)

        for j in range(K):
          
            CI_chosed_val = CI_arr[n == 1]
            CI_chosed_index = np.where(n == 1)[0]

           
            if len(CI_chosed_val) == 0:
                print("No more nodes")
                break

      
            valid_indices = [idx for idx in CI_chosed_index if str(reverse_node_dict[idx]) not in SN_set]

    
            if len(valid_indices) == 0:
                print("No more nodes")
                break


            valid_CI_vals = CI_arr[valid_indices]
            max_index = np.argmax(valid_CI_vals)
            node = valid_indices[max_index]

 
            n[node] = 0
            selected_node = reverse_node_dict[node]


            if str(selected_node) not in SN_set:
                seeds.append(str(selected_node))
            else:
                print(f"Warning")
            b = time.time()
            sumt = 0
            
            self.dataLoader.print_to_file("HCI seeds:",seeds)
            model = HCIC_model(SN, set(seeds), self.vertices, self.hyperedges, alpha, beta)
            for i in range(1000):
                P_number, N_number = model.extend_hyper_IC(vertice_to_hedge_weight)
                sumt += sum(N_number)
            self.dataLoader.print_to_file("The number of HCI seeds:", j, "时,cost time:",b-a,"Negative infected numbers", sumt / 1000)


        return seeds

    
    #celf_greedy
    def celf_greedy(self, K, SN, vertices, hyperedges, vertice_to_hedge_weight, alpha, beta):
        
        SeP_star = set()  
        Queue = []  
        
        I_empty = self.influence_spread(SN, set(), vertices, hyperedges, vertice_to_hedge_weight, alpha, beta)
        self.dataLoader.print_to_file("当SP为空集的时候，阻塞的值为：",I_empty)
        for v in tqdm(vertices - SN):
            
            I_v = self.influence_spread(SN, {v}, vertices, hyperedges, vertice_to_hedge_weight, alpha, beta)
            marginal_gain = I_empty - I_v
            Queue.append((v, marginal_gain))

        
        Queue.sort(key=lambda x: x[1], reverse=True)

        
        va = Queue[0][1]
        SeP_star.add(Queue[0][0])
        Queue.pop(0)
        self.dataLoader.print_to_file("The first node:",SeP_star,"the blocking values:",self.influence_spread(SN, {Queue[0][0]}, vertices, hyperedges, vertice_to_hedge_weight, alpha, beta))

        
        for k in range(2, K + 1):
            check = False
            while not check:
                
                v_max = Queue[0][0]
                SeP = SeP_star.union({v_max})
                I_S = self.influence_spread(SN, SeP, vertices, hyperedges, vertice_to_hedge_weight, alpha, beta)
                
                marginal_gain = I_empty - I_S - va
                Queue[0] = (v_max, marginal_gain)

                
                Queue.sort(key=lambda x: x[1], reverse=True)

                
                if Queue[0][0] == v_max:
                    check = True

            
            va += Queue[0][1]
            SeP_star.add(Queue[0][0])
            self.dataLoader.print_to_file("No",k,",the selected seed set:",SeP_star,"the value:",self.influence_spread(SN, SeP_star, vertices, hyperedges, vertice_to_hedge_weight, alpha, beta))
            Queue.pop(0)

        return SeP_star
    def influence_spread(self, SN, SeP_start, vertices, hyperedges, vertice_to_hedge_weight, alpha, beta):
        
        I_initial = 0
        model = HCIC_model(SN, SeP_start, self.vertices, self.hyperedges, alpha, beta)
        for _ in range(1000):
            
            p_number, n_number = model.extend_hyper_IC(SN, SeP_start, vertices, hyperedges, vertice_to_hedge_weight, alpha, beta)
            I_initial += sum(n_number)  
        I_initial /= 1000  
        return I_initial