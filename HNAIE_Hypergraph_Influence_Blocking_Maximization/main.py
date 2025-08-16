import time

from readfile import DataLoader
from compare_alg import Algs
from hypergraph import Hypergraph_tool
from HCIC_diffusion_model import HCIC_model
from alg_HNAIE import HNAIE





if __name__ == "__main__":
    file_path = "HNAIE\datasets\Restaurants-Rev.txt"

    dataLoader = DataLoader(file_path)
    print("Start import file ",file_path)
    vertices, hyperedges = dataLoader.read_hypergraph_from_txt()
    print(len(vertices))
    print(len(hyperedges))

    tool = Hypergraph_tool(hyperedges, vertices)
    hyperdegrees = tool.compute_hyperdegrees()
    # degrees = tool.compute_degrees()

    weight_dic = {hedge: len(hyperedges[hedge]) for hedge in hyperedges}

    vertice_to_hedge_weight_Hdegree = {v: {hedge: 1/hyperdegrees[v] for hedge in tool.find_hyperedges_containing_node(v)} for v in vertices}
    # vertice_to_hedge_weight_size = {v: {hedge: len(hyperedges[hedge])/sum(tool.cal_temp_weight_key(v, weight_dic).values())
    #                                 for hedge in tool.find_hyperedges_containing_node(v)} for v in vertices}
    
    Algorithm = Algs(hyperedges, vertices, dataLoader)
    node_dict, reverse_node_dict, N, M = Algorithm.generateMap(file_path)
    df_hyper_matrix = Algorithm.changeEdgeToMatrix(file_path)

    #打开存储的文件
    file_in = r"HNAIE\new_results\Restaurant.txt"
    dataLoader.open_output_file(file_in)

    k = 30
    SN = set(tool.get_top_k_nodes_by_degree(k))
    dataLoader.print_to_file(SN)
    alpha = 0.1
    beta = 0.15
    
    #Degree
    dataLoader.print_to_file("Degree method:————————————")
    for d in range(0, 31, 1):
        sumt = 0
        a = time.time()
        SP = set(Algorithm.get_top_k_nodes_excluding_sn(d, SN))
        b = time.time()
        diffusion_model = HCIC_model(SN,SP,vertices,hyperedges,alpha,beta)
        dataLoader.print_to_file("The seeds:",SP)
        
        for i in range(1000):
            P_number, N_number = diffusion_model.extend_hyper_IC(vertice_to_hedge_weight_Hdegree)
            sumt += sum(N_number)
        dataLoader.print_to_file("The number of Degree seeds:", d, "时,cost time:",b-a,"Negative infected numbers", sumt / 1000)
        

    #HDegree
    dataLoader.print_to_file("HDegree method:————————————")
    for d in range(0, 31, 1):
        sumt = 0
        a = time.time()
        SP = set(Algorithm.get_top_k_hyperdegree_nodes_excluding_sn(d, SN))
        b = time.time()
        diffusion_model = HCIC_model(SN,SP,vertices,hyperedges,alpha,beta)
        dataLoader.print_to_file("The seeds:",SP)
        
        for i in range(1000):
            P_number, N_number = diffusion_model.extend_hyper_IC(vertice_to_hedge_weight_Hdegree)
            sumt += sum(N_number)
        dataLoader.print_to_file("The number of HDegree seeds:", d, "时,cost time:",b-a,"Negative infected numbers", sumt / 1000)

    # Random
    dataLoader.print_to_file("Random method:————————————")
    for d in range(0, 31, 1):
        sumt = 0
        a = time.time()
        SP = set(Algorithm.random_select_nodes(d, SN))
        b = time.time()
        diffusion_model = HCIC_model(SN,SP,vertices,hyperedges,alpha,beta)
        dataLoader.print_to_file("The seeds:",SP)
        
        for i in range(1000):
            P_number, N_number = diffusion_model.extend_hyper_IC(vertice_to_hedge_weight_Hdegree)
            sumt += sum(N_number)
        dataLoader.print_to_file("The number of Random seeds:", d, "时,cost time:",b-a,"Negative infected numbers", sumt / 1000)

    #HADP
    dataLoader.print_to_file("HADP method:————————————")
    for d in range(0, 31, 1):
        sumt = 0
        a = time.time()
        SP = set(Algorithm.getSeeds_hdd(df_hyper_matrix, N, d, reverse_node_dict, SN))
        b = time.time()
        diffusion_model = HCIC_model(SN,SP,vertices,hyperedges,alpha,beta)
        dataLoader.print_to_file("The seeds:",SP)
        
        for i in range(1000):
            P_number, N_number = diffusion_model.extend_hyper_IC(vertice_to_hedge_weight_Hdegree)
            sumt += sum(N_number)
        dataLoader.print_to_file("The number of Random seeds:", d, "时,cost time:",b-a,"Negative infected numbers", sumt / 1000)


    #HCI
    dataLoader.print_to_file("HCI_2方法：————————————")
    a = time.time()
    SP = set(Algorithm.getSeeds_ci(2,N,31,df_hyper_matrix,reverse_node_dict,SN,alpha,beta,vertice_to_hedge_weight_Hdegree))
    b = time.time()

    #HRIS
    dataLoader.print_to_file("HRIS method:————————————")
    lamda = 0.1
    num = 200
    for d in range(0, 31, 1):
        sumt = 0
        a = time.time()
        SP = set(Algorithm.getSeeds_ris(N, d, lamda, num, df_hyper_matrix, reverse_node_dict, SN))
        b = time.time()
        diffusion_model = HCIC_model(SN, SP, vertices, hyperedges, alpha, beta)
        dataLoader.print_to_file("The seeds:",SP)
        
        for i in range(1000):
            P_number, N_number = diffusion_model.extend_hyper_IC(vertice_to_hedge_weight_Hdegree)
            sumt += sum(N_number)
        dataLoader.print_to_file("The number of HRIS seeds:", d, "时,cost time:",b-a,"Negative infected numbers", sumt / 1000)

    #HNAIE方法
    dataLoader.print_to_file("HNAIE method:————————————")
    for d in range(0, 31, 1):
        HNAIE_alg = HNAIE(vertices, SN, d, hyperedges, vertice_to_hedge_weight_Hdegree, alpha, beta, weight_dic)
        sumt = 0
        a = time.time()
        SP = set(HNAIE_alg.CMIA_O())
        b = time.time()
        diffusion_model = HCIC_model(SN, SP, vertices, hyperedges, alpha, beta)
        dataLoader.print_to_file("The seeds:",SP)
        
        for i in range(1000):
            P_number, N_number = diffusion_model.extend_hyper_IC(vertice_to_hedge_weight_Hdegree)
            sumt += sum(N_number)
        dataLoader.print_to_file("The number of HNAIE seeds:", d, "时,cost time:",b-a,"Negative infected numbers", sumt / 1000)
    

    #CELF-Greedy方法
    Algorithm.celf_greedy(31, SN, alpha, beta, vertice_to_hedge_weight_Hdegree)


