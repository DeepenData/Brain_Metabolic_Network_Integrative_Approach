def create_toy_metabolism():
    import warnings 
    warnings.filterwarnings("ignore")
    from cobra import Model, Reaction, Metabolite
    model = Model('toy_metabolism')
    # --- Metabolites creation
    A_e    = Metabolite("A[e]")
    A_c    = Metabolite("A[c]")
    B_c    = Metabolite("B[c]")
    atp_c  = Metabolite("atp[c]")
    nadh_c = Metabolite("nadh[c]")
    C_c    = Metabolite("C[c]")
    E_c    = Metabolite("E[c]")
    o2_e   = Metabolite("o2[e]")
    o2_c   = Metabolite("o2[c]")
    E_e    = Metabolite("E[e]")
    # --- Reactions creation
    A_Exchange  = Reaction("A_Exchange")
    A_Uptake    = Reaction("A_Uptake")
    r1          = Reaction("r1")
    r2          = Reaction("r2")
    r3          = Reaction("r3")
    r5          = Reaction("r5")
    OxPhox      = Reaction("OxPhox")
    ATP_demand  = Reaction("ATP_demand")
    C_sink      = Reaction("C_sink")
    O2_Exchange = Reaction("O2_Exchange")
    O2_Uptake   = Reaction("O2_Uptake")
    E_Exchange  = Reaction("E_Exchange")
    E_Uptake    = Reaction("E_Uptake")

    A_Exchange.add_metabolites({A_e: -1})
    A_Uptake.add_metabolites({ A_e: -1, A_c: 1})
    r1.add_metabolites({A_c: -1, atp_c: -1,B_c: 1})
    r2.add_metabolites({B_c: -1, atp_c: 2, nadh_c: .5, C_c: 1})
    r3.add_metabolites({C_c: -1, nadh_c: 4})
    r5.add_metabolites({C_c: -1, nadh_c: -2, E_c: 3})
    OxPhox.add_metabolites({nadh_c: -6, o2_c: -1, atp_c: 2})
    ATP_demand.add_metabolites({atp_c: -1})
    C_sink.add_metabolites({C_c: -1, atp_c: -10})
    O2_Exchange.add_metabolites({o2_e: -1})
    O2_Uptake.add_metabolites({ o2_e: -1, o2_c: 1})
    E_Exchange.add_metabolites({ E_e: -1})
    E_Uptake.add_metabolites({E_e: -1, E_c: 1})

    model.add_reactions([A_Exchange, A_Uptake, r1, r2, r3, r5, OxPhox, \
                        ATP_demand, C_sink, O2_Exchange, O2_Uptake, \
                        E_Exchange, E_Uptake])

    model.reactions.get_by_id("A_Exchange").bounds = (-1000, 0)       
    model.reactions.get_by_id("A_Uptake").bounds = (0, 3.3)     #
    model.reactions.get_by_id("r1").bounds = (0, 1000)          # 
    model.reactions.get_by_id("r2").bounds = (0, 1000)          # 
    model.reactions.get_by_id("r3").bounds = (0, 1000)          # 
    model.reactions.get_by_id("r5").bounds =(-1000, 1000)       # 
    model.reactions.get_by_id("OxPhox").bounds = (0, 1000)      # 
    model.reactions.get_by_id("ATP_demand").bounds = (3, 1000)  # 
    model.reactions.get_by_id("C_sink").bounds = (0, 1000)      # 
    model.reactions.get_by_id("O2_Exchange").bounds = (-1000, 0) # EXC   
    model.reactions.get_by_id("O2_Uptake").bounds = (0, 10)     # 
    model.reactions.get_by_id("E_Exchange").bounds = (-10000, 1000) # EXC     
    model.reactions.get_by_id("E_Uptake").bounds = (-1000, 0) #

    model.reactions.get_by_id("r5").bounds = (0, 1000)
    #model.reactions.get_by_id("E_Exchange").bounds = (0, 1000)
    model.reactions.get_by_id("E_Uptake").bounds = (-1000, 0)
    model.objective = "C_sink"
    return model

def str_match(node_list, pattern):
    import re
    hh  = [re.findall(pattern, a_node) for a_node in node_list]
    nodes = [item for sublist in hh for item in sublist]
    return nodes

def generate_rbd_graph(A, B, p):
    from networkx.algorithms.bipartite.generators import random_graph
    import networkx as nx
    
    import random as rd
    import numpy as np
    rbd_graph = random_graph(A,B,p,directed=False)
    rbd_graph = nx.DiGraph(rbd_graph).to_directed()
    ebunch    = rd.sample(list(rbd_graph.edges),B)
    rbd_graph.remove_edges_from(ebunch)

    while np.invert(nx.is_connected(nx.Graph(rbd_graph))):
        rbd_graph = random_graph(A,B,p,directed=False)
        rbd_graph = nx.DiGraph(rbd_graph).to_directed()
        ebunch    = rd.sample(list(rbd_graph.edges),B)
        rbd_graph.remove_edges_from(ebunch)
    
    print("Bipartite:",
    nx.is_bipartite(rbd_graph),"\nDirected:",
    nx.is_directed(rbd_graph),"\nConnected:",
    nx.is_connected(nx.Graph(rbd_graph)))

    return rbd_graph


def plot_bipartite(G, node_color = "grey"):
    from networkx import bipartite
    import networkx as nx
    import matplotlib.pyplot as plt

    X, Y = bipartite.sets(G)
    pos = dict()
    pos.update( (n, (1, i)) for i, n in enumerate(X) ) # put nodes from X at x=1
    pos.update( (n, (2, i)) for i, n in enumerate(Y) ) # put nodes from Y at x=2
    nx.draw(G,  node_color= node_color,  pos=pos, cmap=plt.get_cmap('viridis'), with_labels=True, font_color='red')
    return plt.show()

def set_node_names(G):
    import networkx as nx
    from networkx import bipartite

    A, B = bipartite.sets(G)

    rxns      = ["rxn" + str(i) for i in range(0,len(A))]
    mets      = ["met" + str(i) for i in range(0,len(B))]
    mapping_A = dict(zip(A, rxns ))
    mapping_B = dict(zip(B,  mets ))

    G = nx.relabel_nodes(G, mapping_A)
    G = nx.relabel_nodes(G, mapping_B)
    return G

def make_union(G1, H1):
    import networkx as nx
    from networkx.algorithms.operators.binary import union
    G1_H1 = union(G1,H1, rename=('G-','H-') )
    G1_H1 = nx.DiGraph(G1_H1).to_directed()
    return G1_H1




def make_bridges(G1_H1):
    degrees   = dict(list(G1_H1.degree))
    node_list = list({key:value for (key, value) in degrees.items() if value <= 2}.keys())      
    
    import networkx as nx 

    def get_bridges(node_list = node_list):   

        G_rxn = str_match(node_list, "G.rxn.*")
        G_met = str_match(node_list, "G.met.*")
        H_rxn = str_match(node_list, "H.rxn.*")
        H_met = str_match(node_list, "H.met.*")

        import random as rd

        G_rxn_samples = rd.sample(G_rxn,2)
        G_met_samples = rd.sample(G_met,2)

        H_rxn_samples = rd.sample(H_rxn,2)
        H_met_samples = rd.sample(H_met,2)

        bridges = list()
        for G_r, G_m, H_r, H_m in zip(G_rxn_samples, G_met_samples, H_rxn_samples, H_met_samples):
                
                e1 = (G_r, H_m,{})
                e2 = (H_m, G_r,{})
                e3 = (H_r, G_m,{})
                e4 = (G_m, H_r,{})
                bridges.append([e1, e2, e3, e4])
                
        bridges = [item for sublist in bridges for item in sublist]
        return bridges
    
    bridges = get_bridges()
    G1_H1.add_edges_from(bridges)  

    import networkx as nx

    print("Bipartite:",
    nx.is_bipartite(G1_H1),"\nDirected:",
    nx.is_directed(G1_H1),"\nConnected:",
    nx.is_connected(nx.Graph(G1_H1)))
    return G1_H1

####################33
###############
from retrying import retry
@retry(stop_max_delay=60000, wait_fixed=0)
def get_two_core_directed_random(A = 16, B = 12, p = 0.11):
    G1    =  generate_rbd_graph(A, B , p ) #A y B, son la cantidad de nodos en cada partición. p: probabilidad de edge
    #Se asignan nombres de pseudo-reacciones y metabolitos (set_node_names).
    G1    =  set_node_names(G1)
    H1    =  generate_rbd_graph(A , B , p)
    H1    =  set_node_names(H1)
    #La union (make_union) de de los grafos y conexión (make_bridges) de los dos grafos (G1 y H1), 
    # cada nombre de nodo tendrá un prefijo de su grafo de origen (G- o H-). 
    # Los grafos de origen son las subredes.
    G1_H1 =  make_union(G1, H1)
    G1_H1 =  make_bridges(G1_H1)
    return G1, G1_H1


from graph_tool.all import *

def set_subgraphs_attr_bipartite(gpath, figure_path):

    import networkx as nx

    NX = nx.read_graphml(gpath)
    GT = load_graph(gpath)
    blocks = minimize_blockmodel_dl(GT).get_blocks()
    graph_draw(GT, vertex_fill_color = "lightblue", output= figure_path,vertex_shape = GT.vp["bipartite"],\
        vertex_size=20, bg_color = "white", edge_pen_width=2.5, edge_color = "grey", edge_marker_size = 10)
    subgra_attrs = list(blocks.a)
    NX = list2attr(NX, list(NX.nodes), "subgraph", subgra_attrs)
    nx.write_graphml(NX, gpath)
    return NX
    

def draw_bipartite(graphml_path, pdf_file ):
    
    G2_H2 = load_graph(graphml_path)
    partitions  = list(G2_H2.vp["bipartite"].a)

    node_names = [G2_H2.vp._graphml_vertex_id[i] for i in range(len(G2_H2.get_vertices()))]
    subgraph = str_match(node_names, "G|H")
    subgraph1 = [ord(subgraph) for subgraph in subgraph]
    vprop = G2_H2.new_vertex_property("int",subgraph1)
    G2_H2.vp["subgraph"] = vprop
    graph_draw(G2_H2,  vertex_shape = G2_H2.vp["bipartite"],\
            vertex_size=15, edge_pen_width=2,vertex_fill_color = G2_H2.vp["subgraph"], \
                output=pdf_file)
from graph_tool.all import *


def set_subgraphs_attr(gpath, figure_path):

    import networkx as nx

    NX = nx.read_graphml(gpath)
    GT = load_graph(gpath)
    blocks = minimize_blockmodel_dl(GT).get_blocks()
    graph_draw(GT, vertex_fill_color = blocks, output= figure_path,  vertex_size=30, bg_color = "white", edge_pen_width=1.5, edge_color = "grey", edge_marker_size = 10)
    subgra_attrs = list(blocks.a)
    NX = list2attr(NX, list(NX.nodes), "subgraph", subgra_attrs)
    nx.write_graphml(NX, gpath)
    return NX


def draw_projection(graphml_path, pdf_file ):
    
    G2_projected = load_graph( graphml_path)
    node_names = [G2_projected.vp._graphml_vertex_id[i] for i in range(len(G2_projected.get_vertices()))]
    subgraph = str_match(node_names, "G|H")
    subgraph1 = [ord(subgraph) for subgraph in subgraph]
    vprop = G2_projected.new_vertex_property("int",subgraph1)
    G2_projected.vp["subgraph"] = vprop
    graph_draw(G2_projected, \
            vertex_size=15, edge_pen_width=2,vertex_fill_color = G2_projected.vp["subgraph"], \
                output=pdf_file)

def make_rxn_projection(G):
    from networkx import bipartite
    from networkx import projected_graph
    import networkx as nx
    rxns, mets = bipartite.sets(G)
    projected = projected_graph(nx.Graph(G), rxns)
    return projected

from retrying import retry
@retry(stop_max_delay=60000, wait_fixed=0)
def get_two_core_directed_FreeScale(A, p):
    G2    =  generate_pa_graph(A, p)
    G2    =  set_node_names(G2)
    H2    =  generate_pa_graph(A , p )
    H2    =  set_node_names(H2)
    G2_H2 =  make_union(G2, H2)
    G2_H2 =  make_bridges(G2_H2)
    return G2, G2_H2 


def plot_bipartite_subgraphs(G1_H1):
    from networkx import bipartite
    import networkx as nx
    import matplotlib.pyplot as plt

    X, Y = bipartite.sets(G1_H1)

    GH=str_match(list(G1_H1.nodes), "G|H") 

    import re

    GH = [re.sub('G','green', gh) for gh in GH ]
    GH = [re.sub('H','blue', gh) for gh in GH ]

    pos = dict()
    pos.update( (n, (1, i)) for i, n in enumerate(X) ) # put nodes from X at x=1
    pos.update( (n, (2, i)) for i, n in enumerate(Y) ) # put nodes from Y at x=2
    nx.draw_networkx(G1_H1, arrowsize = 5, width = .6, node_size =30, alpha = 0.6, \
        node_color= GH,  pos=pos, edge_color = "grey",  \
        with_labels=False, font_color='red')
    return plt.show()

def generate_pa_graph(A, p):
    from networkx.algorithms.bipartite.generators import preferential_attachment_graph
    import networkx as nx
    import random as rd
    import numpy as np
    from networkx import bipartite
    pa_graph = preferential_attachment_graph(aseq = A,  p = p)
    pa_graph = nx.DiGraph(pa_graph).to_directed()
    ebunch    = rd.sample(list(pa_graph.edges),len(A))
    pa_graph.remove_edges_from(ebunch)
    

    while  np.invert(nx.is_connected(nx.Graph(pa_graph))) :
        pa_graph = preferential_attachment_graph(aseq = A,  p = p)
        pa_graph = nx.DiGraph(pa_graph).to_directed()
        ebunch    = rd.sample(list(pa_graph.edges),len(A)-2)
        pa_graph.remove_edges_from(ebunch)
        #A2, B2 = bipartite.sets(pa_graph)
    A2, B2 = bipartite.sets(pa_graph)
    print("Bipartite:",
    nx.is_bipartite(pa_graph),"\nDirected:",
    nx.is_directed(pa_graph),"\nConnected:",
    nx.is_connected(nx.Graph(pa_graph)), "\nlen(A2) > len(B2)" ,(len(A2) > len(B2)))

    return pa_graph

def plot_degree_distribution(G):
    import collections
    import matplotlib.pyplot as plt
    import networkx as nx

    G = nx.Graph(G)

    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
    degreeCount = collections.Counter(degree_sequence)
    deg, cnt = zip(*degreeCount.items())

    fig, ax = plt.subplots()
    plt.bar(deg, cnt, width=0.80, color="b")

    plt.title("Degree Histogram")
    plt.ylabel("Count")
    plt.xlabel("Degree")
    ax.set_xticks([d + 0.4 for d in deg])
    ax.set_xticklabels(deg)

    # draw graph in inset
    plt.axes([0.4, 0.4, 0.5, 0.5])
    Gcc = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
    pos = nx.spring_layout(G)
    plt.axis("off")
    nx.draw_networkx_nodes(G, pos, node_size=20)
    nx.draw_networkx_edges(G, pos, alpha=0.4)
    return plt.show()

def get_largest_component(grafo): 
    import networkx as nx
    largest_component = max(nx.connected_components(grafo), key=len)
    G = grafo.subgraph(largest_component)
    return G

def cobra_to_bipartite(model):
    import networkx as nx
    from   cobra.util.array import create_stoichiometric_matrix
    import numpy as np
    from sklearn.preprocessing import Binarizer
    import warnings 
    from scipy.sparse import csr_matrix
    from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix


    warnings.filterwarnings("ignore")
    #extraer matriz estequiométrica
    S_matrix = create_stoichiometric_matrix(model)
    #convertir todas las entradas valores positivos
    S_matrix = abs(S_matrix)
    #binarizar
    S_matrix = Binarizer().fit_transform(S_matrix)

    S_matrix = csr_matrix(S_matrix) # Convierte la matriz a una matriz dispersa
    grafo = from_biadjacency_matrix(S_matrix) # Usa la dispersa para un grafo bipartito
    
    metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    metabolites_n     = len(model.metabolites) # Numero de metabolitos
    metabolites_names =    [model.metabolites[i].id for i in range(0, metabolites_n) ]

    reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones
    reactions_n       = len(model.reactions)   # Numero de reacciones
    reactions_names   =    [model.reactions[i].id   for i in range(0, reactions_n)   ]

    names_mapped =  dict(zip( metabolites_nodes + reactions_nodes, metabolites_names + reactions_names))
    grafo = nx.relabel_nodes(grafo, names_mapped)
   # grafo = get_largest_component(grafo)
    return grafo


def cobra_to_networkx_rxn_projection(modelo):
    import networkx as nx
    from   cobra.util.array import create_stoichiometric_matrix
    import numpy as np
    from sklearn.preprocessing import Binarizer
    import warnings 
    warnings.filterwarnings("ignore")

    assert str(type(modelo)) == "<class 'cobra.core.model.Model'>", "El objeto debe ser un modelo, no un optimizado (modelo.optimize())"
    #extraer matriz estequiométrica
    S_matrix = create_stoichiometric_matrix(modelo)
    #convertir todas las entradas valores positivos
    S_matrix = (abs(S_matrix) )
    #transformar a enteros 
    S_matrix = S_matrix.astype(np.int)
    #Multiplicacion por la derecha para proyectar en el espacio de las reacciones
    projected_S_matrix = np.matmul(S_matrix.T, S_matrix)
    #rellenar diagonal con ceros
    np.fill_diagonal(projected_S_matrix, 0) 
    #binarizar
    projected_S_matrix = Binarizer().fit_transform(projected_S_matrix)
    #crear grafo networkx
    G = nx.convert_matrix.from_numpy_matrix( projected_S_matrix )
    #hacer diccionario con los nombres de las reacciones
    node_dict   = lambda l : dict( zip( list(G.nodes), l ) )
    reaction_dict = node_dict( [reaction.id for reaction in modelo.reactions] )
    #Renombrar los nodos usando el diccionario
    G = nx.relabel_nodes(G, reaction_dict, copy=True) # Revisar que esto este antes de la remoción del grafo
    #G = get_largest_component(G)
    return G


def list2attr(grafo, nodos, nombre, atributos):
    """Toma dos listas: nombres de nodos y atributos; y las añade a un grafo

    Parameters
    ----------
    grafo : bipartite_graph
        Un grafo bipartito de NetworkX 
    nodos: list
        Una lista de nombres de nodos
    nombre : str
        Nombre del atributo siendo asignado. Ie. nombre de la columna en Gephi
    atributos : list
        Una lista de valores para los nodos en la lista anterior.
        Ambas listas deben ser del mismo largo. 
    
    Returns
    -------
    grafo: bipartite_graph
        Un grafo bipartito con un nuevo atributo para un set de nodos "nodos". 
    """
    assert len(nodos) == len(atributos), "Ambas listas deben ser del mismo largo."
    tmp_list = { nodos[i] : atributos[i] for i in range(0, len(atributos)) }
    
    from networkx import set_node_attributes
    set_node_attributes(grafo, tmp_list, nombre ) # Añade los atributos
    
    return grafo

def solution2attr(solucion, grafo, estandarizar=False, umbral=1e-7):
    """Docstring

    Parameters
    ----------
    solucion : 
    grafo : 
    estandarizar : bool, default False
        Define si se aplican estandarizaciones al modelo final. Estas son entre
    umbral : float
        El umbral en que un flujo o sensibilidad se consideran 0. 
    """

    metabolites_nodes = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 0] # Crea una lista de metabolitos
    reactions_nodes   = [n for n, d in grafo.nodes(data=True) if d["bipartite"] == 1] # Crea una lista de reacciones

    print("Metabolitos:", len(metabolites_nodes), " | Reacciones:", len(reactions_nodes)) # Debug

    from numpy import abs # Permite absolutos vectorizados

    shadow_prices = abs(solucion.shadow_prices).tolist() # De solucion, pasa a lista
    fluxes        = abs(solucion.fluxes).tolist()       # posiblemente es más rapido de usar que desde el modelo
    reduced_costs = abs(solucion.reduced_costs).tolist() # considerando el menor parsing pero si requiere pandas

    print("Shadow prices:", len(shadow_prices),"| Fluxes:", len(fluxes),"| Reduced costs:", len(reduced_costs)) # Debug

    # Neutraliza sub-umbral (por defecto 0.0000001)
    shadow_prices = [ 0 if i < umbral else i for i in shadow_prices]
    fluxes        = [ 0 if i < umbral else i for i in fluxes]
    reduced_costs = [ 0 if i < umbral else i for i in reduced_costs]

    if estandarizar == True:
        def estandariza(tmp):
            from numpy import array, log10, inf, min, max
            tmp = log10(array(tmp))
            tmp[tmp == -inf] = tmp[tmp != -inf].min()
            tmp = (tmp - tmp.min() ) / ( tmp.max() - tmp.min() )
            return tmp
        
        shadow_prices = estandariza(shadow_prices)
        fluxes = estandariza(fluxes)
        reduced_costs = estandariza(reduced_costs)

    from networkx import set_node_attributes

    # ASIGNA Shadow_prices, fluxes, reduced_costs
    set_node_attributes(grafo, { metabolites_nodes[i] : shadow_prices[i] for i in range(0, len(shadow_prices)) } , "Shadow Price") 
    set_node_attributes(grafo, { reactions_nodes[i] : fluxes[i] for i in range(0, len(fluxes)) } , "Flux") 
    set_node_attributes(grafo, { reactions_nodes[i] : reduced_costs[i] for i in range(0, len(reduced_costs)) } , "Reduced cost") 

    #grafo.nodes(data=True) # Debug

    return grafo
