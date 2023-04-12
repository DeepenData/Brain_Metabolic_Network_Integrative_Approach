# %% --- 
import sys

from numpy.core.numeric import ones
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
supp_path    = local_working_directory + "/Supplementary/"
supp_path_results    = local_working_directory + "/Supplementary/" + "results/"
from aa_functions import *
import numpy as np
from cobra import *
A_uptake_range  = np.linspace(0, 20, num=50)
O2_uptake_range = np.linspace(0, 10, num=50)
model = []
model = create_toy_metabolism()

obj = np.empty([len(A_uptake_range), len(O2_uptake_range)])
A = np.empty([len(A_uptake_range), len(O2_uptake_range)])
O2 =np.empty([len(A_uptake_range), len(O2_uptake_range)])

for i, A_Uptake in enumerate(A_uptake_range):
  for j, O2_Uptake in enumerate(O2_uptake_range):
    model.reactions.get_by_id("A_Uptake").bounds = (A_Uptake, A_Uptake)
    model.reactions.get_by_id("O2_Uptake").bounds = (O2_Uptake, O2_Uptake)

    sol       = model.optimize()

    if sol.objective_value <= 0:
      obj[i,j]  = 0
    else:
      obj[i,j]  = sol['C_sink']
      A[i,j]  = sol['A_Uptake']
      O2[i,j] = sol['O2_Uptake']
X = A
Y = O2
Z = obj
import pandas as pd

s = abs(pd.read_csv(supp_path_results + 'samples_leftraru.csv', header=None).T)
cols_s = ['A_Exchange','A_Uptake','r1','r2','r3','r5','OxPhox','ATP_demand',\
    'C_sink','O2_Exchange','O2_Uptake','E_Exchange','E_Uptake']
s.columns = list(cols_s)

sX   = s.A_Uptake
sY   = s.O2_Uptake
sZ   = s.C_sink

import plotly.graph_objects as go
fig = go.Figure(data=[ go.Surface(x=X, y= Y, z=Z, opacity=0.8),\
  go.Scatter3d(x=sX, y= sY, z=sZ, mode ='markers', 
                                   marker = dict(
                                     size = 3,
                                     colorscale ='Viridis',
                                     opacity = 0.5
                                   ))])
fig.update_layout(coloraxis_showscale=False, scene=dict(
        xaxis_title='A_Uptake',
        yaxis_title='O2_Uptake',
        zaxis_title='C_sink'),margin=dict(l=0, r=100, b=0, t=0)
)                                   
fig.show()
fig.write_image(supp_path_results +"samples_and_phpp.png")
# %% --- 
#import IPython
#IPython.Application.instance().kernel.do_shutdown(True) 
import numpy as np
import pandas as pd
A_uptake_range  = np.linspace(0, 20, num=50)
O2_uptake_range = np.linspace(0, 10, num=50)
model.reactions.get_by_id("A_Uptake").bounds = ( min(A_uptake_range), max(A_uptake_range))
model.reactions.get_by_id("O2_Uptake").bounds = (min(O2_uptake_range), max(O2_uptake_range))

model.reactions.get_by_id("A_Exchange").bounds = ( -max(A_uptake_range), max(A_uptake_range))
model.reactions.get_by_id("O2_Exchange").bounds = (-max(O2_uptake_range), max(O2_uptake_range))

sol = model.optimize()

fba0 = pd.DataFrame([( r.reaction, r.id, r.bounds) for r in model.reactions], list(sol.fluxes))
fba = pd.DataFrame([( r.reaction, r.id, r.bounds) for r in model.reactions])
fba.columns = ["Reaction","Name","Bounds"]
fba.reset_index(drop=True, inplace=True)

fba.index = list(range(1,14))
import dataframe_image as dfi
dfi.export(fba, supp_path_results + 'model_rxns.png')
# %% --- 
import IPython
IPython.Application.instance().kernel.do_shutdown(True) 
import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
supp_path    = local_working_directory + "/Supplementary/"
supp_path_results    = local_working_directory + "/Supplementary/" + "results/"

import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
supp_path    = local_working_directory + "/Supplementary/"
from aa_functions import *
import numpy as np
from cobra import *
A_uptake_range  = np.linspace(0, 20, num=50)
O2_uptake_range = np.linspace(0, 10, num=50)
model = []
model = create_toy_metabolism()

obj = np.empty([len(A_uptake_range), len(O2_uptake_range)])
A = np.empty([len(A_uptake_range), len(O2_uptake_range)])
O2 =np.empty([len(A_uptake_range), len(O2_uptake_range)])

for i, A_Uptake in enumerate(A_uptake_range):
  for j, O2_Uptake in enumerate(O2_uptake_range):
    model.reactions.get_by_id("A_Uptake").bounds = (A_Uptake, A_Uptake)
    model.reactions.get_by_id("O2_Uptake").bounds = (O2_Uptake, O2_Uptake)

    sol       = model.optimize()

    if sol.objective_value <= 0:
      obj[i,j]  = 0
    else:
      obj[i,j]  = sol['C_sink']
      A[i,j]  = sol['A_Uptake']
      O2[i,j] = sol['O2_Uptake']
X = A
Y = O2
Z = obj

model.reactions.get_by_id("A_Uptake").bounds = (0, 10.4)
model.reactions.get_by_id("O2_Uptake").bounds = (0, 1000)

sol       = model.optimize()
C_sink = sol.to_frame().loc['C_sink',:].fluxes
A_Uptake = sol.to_frame().loc['A_Uptake',:].fluxes
O2_Uptake = sol.to_frame().loc['O2_Uptake',:].fluxes

import plotly.graph_objects as go
import numpy as np
fig = go.Figure(data=[ 
      go.Surface(x=X, y= Y, z=Z, opacity=0.4,
          contours = {
            "z": {"show": True}
    }),\
      go.Scatter3d(x=  np.array(A_Uptake), 
               y=  np.array(O2_Uptake), 
               z=  np.array(C_sink), mode ='markers', 
                                   marker = dict(
                                     size = 10,
                                     #color = df['petal_width'],
                                     colorscale ='Viridis',
                                     opacity = 1
                                   ))])

fig.update_layout(coloraxis_showscale=False, scene=dict(
        xaxis_title='A_Uptake',
        yaxis_title='O2_Uptake',
        zaxis_title='C_sink',
    ),margin=dict(l=0, r=100, b=0, t=0))

fig.show()   
fig.write_image(supp_path_results +"phpp_optimal_solution.png")
# %% --- 
# 
#import IPython
#IPython.Application.instance().kernel.do_shutdown(True) 
import sys
local_working_directory = '/home/alejandro/AA_PostDoc/Acevedo_et_al_2021_neuron_astrocyte_metabolism'
sys.path.insert(0,local_working_directory)
supp_path    = local_working_directory + "/Supplementary/"
supp_path_results    = local_working_directory + "/Supplementary/" + "results/"
from aa_functions import *
import numpy as np
from cobra import *
import networkx as nx
from   cobra.util.array import create_stoichiometric_matrix
from scipy.sparse import csr_matrix
from networkx.algorithms.bipartite.matrix import from_biadjacency_matrix
import pandas as pd
from graph_tool import *
from cobra import Model, Reaction, Metabolite
import sys
figure1_path = local_working_directory + "/Figure_1"
figure2_path = local_working_directory + "/Figure_2"
figure3_path = local_working_directory + "/Figure_3"
figure4_path = local_working_directory + "/Figure_4"

model = create_toy_metabolism()
S_matrix = create_stoichiometric_matrix(model)
S_matrix = csr_matrix(S_matrix) # Convierte la matriz a una matriz dispersa
nx_biSmallMetab = from_biadjacency_matrix(S_matrix, create_using = nx.DiGraph) # Usa la dispersa para un grafo bipartito
nx_biSmallMetab_2 = nx.DiGraph((cobra_to_bipartite(model)))
edges_to_remove = list(nx_biSmallMetab_2.edges)
nx_biSmallMetab_2.remove_edges_from(edges_to_remove)
S_matrix = create_stoichiometric_matrix(model)
S_matrix = csr_matrix(S_matrix)
r,c =S_matrix.nonzero()

met_ids = [m.id for m in model.metabolites]
mets = [met_ids[row] for row in r]

rxn_ids = [r.id for r in model.reactions]
rxns = [rxn_ids[col] for col in c]

coef = list(np.array(S_matrix[r,c]).flatten())

my_edges = list()
for m, r, c in zip(mets, rxns, coef):
  if c < 0:
    my_edges.append((m,r))
  else:
    my_edges.append((r,m))


my_edges.append(('O2_Exchange', 'o2[e]'))
my_edges.append(('A_Exchange', 'A[e]'))
my_edges.append(('E_Exchange', 'E[e]'))


nx_biSmallMetab_2.add_edges_from(my_edges)
nx_biSmallMetab_2.edges(data=True)

shapes = \
dict(zip(list(nx_biSmallMetab_2.nodes),1*(np.array(list(nx.get_node_attributes(nx_biSmallMetab_2, "bipartite").values()))==0) ))
nx.set_node_attributes(nx_biSmallMetab_2, shapes, "shapes")


from networkx import bipartite

nx_biSmallMetab_3 = nx_biSmallMetab_2
mets, rxns  = bipartite.sets(nx_biSmallMetab_3)

import numpy as np


#A_uptake_range  = np.linspace(0, 20, num=50)
#O2_uptake_range = np.linspace(0, 10, num=50)
#mid = round(len(A_uptake_range)/2)
#model.reactions.get_by_id("A_Uptake").bounds =  0, A_uptake_range[mid] 
model.reactions.get_by_id("A_Uptake").bounds = (0, 10.4)
model.reactions.get_by_id("O2_Uptake").bounds = (0, 1000)

model.objective = "C_sink"
sol       = model.optimize()

fba = pd.DataFrame([( r.reaction, r.id, r.bounds) for r in model.reactions], list(sol.fluxes))

fluxes = \
dict(zip( list(fba.iloc[:,1])  , round(sol.fluxes,2).astype("str")  ))


fluxes_mets = dict(zip(mets, [''] * 10))

fluxes.update(fluxes_mets)

nx.set_node_attributes(nx_biSmallMetab_3, fluxes, "fluxes")

rd = list(abs(round(sol.reduced_costs,2)).astype("str") )
fl = list(round(sol.fluxes,2).astype("str"))
rxn_num = list(np.linspace(1,len(fl),len(fl)).astype("int").astype("str"))
hola = list()
for n, x, y in zip(rxn_num, fl,  rd):
    hola.append(n + ": " + x + ", "+ y)

sensi = \
dict(zip( list(fba.iloc[:,1])  , hola ))
sensi_mets = dict(zip(mets, mets))#[''] * 10))

sensi.update(sensi_mets)

nx.set_node_attributes(nx_biSmallMetab_3, sensi, "sensi")

nx_biSmallMetab_3.nodes(data=True)

nx_biSmallMetab_3 = list2attr(nx_biSmallMetab_3,list(nx_biSmallMetab_3.nodes()), "position", list(np.ones(23).astype("int")) )

nx_biSmallMetab_3.nodes(data=True)
# %% --- 

nx.write_graphml(nx_biSmallMetab_3, supp_path_results + "Bipartite_small_metabolism.graphml")

#nx_biSmallMetab_2.nodes(data=True)
#plot_bipartite(nx_biSmallMetab_2)
gpath = supp_path_results +"Bipartite_small_metabolism.graphml"
GT    = load_graph(gpath)
#GT.list_properties()
figure_path = supp_path_results +"Bipartite_small_metabolism.png"
import numpy as np              
import matplotlib

pos = sfdp_layout(GT)


graph_draw(GT, pos=pos,vertex_font_size = 10,vcmap=matplotlib.cm.tab20, output_size=(1400, 1400), \
  vertex_fill_color = GT.vertex_properties['shapes'], #vertex_text_position =  GT.vertex_properties['position'],\
   output= figure_path,vertex_shape = GT.vertex_properties['shapes'],adjust_aspect = True, fit_view_ink = True,\
                       vertex_text = GT.vertex_properties['sensi'], \
           vertex_size= 60, bg_color = "white", edge_pen_width=5.5, edge_color = "grey", edge_marker_size = 12)                           