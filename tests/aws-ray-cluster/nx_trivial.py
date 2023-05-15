"""
Este script simplemente lanza actores distribuidos para procesar un grafo grande
de networkx. La idea de usar actores es que permite crear objetos state-full, 
versus solo usar funciones state-less, como en las primeras iteraciones de Ray. 
"""

import ray
import random
import networkx as nx
import time 

@ray.remote(num_cpus=0.25)
class RandomGraphActor:

    # Constructor que acepta un grafo    
    def __init__(self, G):
        self.G = G

    # Metodo para computar centralidad entre dos nodos
    def get_connectivity(self, nodo_s, nodo_t):
        time.sleep(2.5) # Procesamiento lento
        return nx.node_connectivity(self.G, nodo_s, nodo_t)

     # Los actores solo permiten usar metodos para preguntar informacion interna
    def get(self):
        return self.G

    def n(self):
        return nx.number_of_nodes( self.G )
    

# Construye un grafo grande y lo deja en el almacen de objetos de Ray
# para asi evitar desperdiciar IO con cada instantacion de los actores
grafo_absurdo = nx.fast_gnp_random_graph(1000, 0.25) # grafo ...
g_absurdo_remoto = ray.put( grafo_absurdo )         # ... a object_store


print("CREANDO ACTORES")

num_actors = 50

actors = [] # <- Lista de actores
for _ in range(num_actors):
    # n = random.randint(min_nodes, max_nodes)
    actor = RandomGraphActor.remote(g_absurdo_remoto)
    actors.append(actor)

print("\nLANZANDO PROCESAMIENTOS\n")

min_nodes = 100
max_nodes = 1000

futures = [] 
for actor in actors:
    s, t = random.randint(min_nodes, max_nodes), random.randint(min_nodes, max_nodes)
    print(f"Graph has { ray.get(actor.n.remote()) } nodes. Path {s} to {t}")
    
    # Lanza una lista de tareas de computo para cada actor, llamando al metodo
    # que traza centeralidad entre nodos. Esto es state-full. 
    futures.append( 
        actor.get_connectivity.remote(s,t) 
    )
    # En este momento, la lista existe como una tarea asincrona en el clsuter

print("\DOING WHATEVER ELSE\n")
ray.cluster_resources()

print("\nRESULTADOS\n")
for future in futures:
    print(f"Connectivity of the graph: {ray.get(future)}")

# TODO: revivir codigo de Namu para nombres ordenados de tareas, el uso de while,. 
#   y los dict-comprehension con nombres y argumentos de funciones. 