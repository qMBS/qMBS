#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import time
from qiskit import *
from qiskit import QuantumCircuit, transpile
from qiskit.compiler import assemble
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
from qiskit import IBMQ
from qiskit import Aer
import pandas as pd
import networkx as nx
import random

def downsample_bipartite_graph(data, target_min=10, target_max=20):
    G = nx.Graph()
    G.add_edges_from(data.values)

    while len(G.nodes) > target_max:
        degrees = dict(G.degree())
        total_degree = sum(degrees.values())
        if total_degree == 0:
            break

        probabilities = {node: degree / total_degree for node, degree in degrees.items()}
        sampled_nodes = set()

        for node in G.nodes:
            if random.random() < probabilities[node]:
                sampled_nodes.add(node)
                sampled_nodes.update(G.neighbors(node))

        if len(sampled_nodes) < target_min:
            sampled_nodes.update(random.sample(list(G.nodes), min(target_min, len(G.nodes))))

        G = G.subgraph(sampled_nodes).copy()

        if len(G.nodes) <= target_min:
            break

    return G




def bipartite_to_adjacency_matrix(edges, v_nodes, u_nodes):
    v_index = {v: i for i, v in enumerate(v_nodes)}
    u_index = {u: i for i, u in enumerate(u_nodes)}
    adjacency_matrix = np.zeros((len(v_nodes), len(u_nodes)), dtype=int)
    for v, u in edges:
        if v in v_index and u in u_index:
            adjacency_matrix[v_index[v], u_index[u]] = 1
    return adjacency_matrix

def generate_bipartite_graph_variables(n, m):
    for i in range(1, n + 1):
        exec(f"global v{i}; v{i} = {i - 1}")
    for j in range(1, m + 1):
        exec(f"global u{j}; u{j} = {n + j - 1}")
    total_edges = n * m
    for k in range(1, total_edges + 1):
        exec(f"global e{k}; e{k} = {n + m + k - 1}")
        exec(f"global e{k}_; e{k}_ = {n + m + total_edges + k - 1}")
    exec(f"global bic; bic = {n + m + total_edges * 2}")
    cl_index = n + m + total_edges * 2 + 1
    for i in range(1, n + 1):
        for j in range(0, i + 1):
            exec(f"global cl{i}{j}; cl{i}{j} = {cl_index}")
            cl_index += 1
    cr_index = cl_index
    for i in range(1, n + 1):
        for j in range(0, i + 1):
            exec(f"global cr{i}{j}; cr{i}{j} = {cr_index}")
            cr_index += 1
    exec(f"global ce1; ce1 = {cr_index}")
    exec(f"global ce2; ce2 = {cr_index + 1}")
    exec(f"global ce4; ce4 = {cr_index + 2}")
    exec(f"global o; o = {cr_index + 3}")

################ downsample ##################
# file_path = '/Users/xl/Downloads/youtube-groupmemberships.txt'
# #download from http://konect.cc/networks/

# data = pd.read_csv(file_path, sep='\t', header=None, names=["V", "U"])
# data_cleaned = data[data["V"].str.match(r"^\d+ .*", na=False)]
# data_cleaned = data_cleaned["V"].str.split(expand=True).rename(columns={0: "V", 1: "U"})
# data_cleaned["V"] = data_cleaned["V"].astype(int)
# data_cleaned["U"] = data_cleaned["U"].astype(int)

# downsampled_graph = downsample_bipartite_graph(data_cleaned, target_min=10, target_max=20)
# downsampled_edges = list(downsampled_graph.edges())

# output_file_path = '/Users/xl/Downloads/down_youtube-groupmemberships.txt'
# with open(output_file_path, 'w') as f:
#     for edge in downsampled_edges:
#         f.write(f"{edge[0]}\t{edge[1]}\n")
################ downsample ##################

D = np.array([[1,0],[1,1],], dtype=bool)
n, m = 2, 2
generate_bipartite_graph_variables(n, m)

virtues = [globals()[f"e{i}_"] for i in range(1, n * m + 1)]
# APIkey = "pFb_6Nqdn80VEutd7pluKBzrySJLXZi"

# The user needs to apply for their own token on the IBM platform.
TOKEN = "********"



def Statepre(*qubits, QuantumCircuit):
    circuit = QuantumCircuit
    for q in qubits[:-1]:  
        circuit.h(q)
    circuit.x(qubits[-1])  
    circuit.h(qubits[-1])

    
def C_CX(ho,so,targ,QuantumCircuit):
    QuantumCircuit.x(ho)
    QuantumCircuit.ccx(ho,so,targ)
    QuantumCircuit.x(ho)



def CCCZ(c1, c2, c3, targ, QuantumCircuit):
    circuit = QuantumCircuit
    ce1 = circuit.num_qubits - 3
    ce2 = circuit.num_qubits - 2
    circuit.ccx(c1, c2, ce1)
    circuit.ccx(c3, ce1, ce2)
    circuit.cz(targ, ce2)
    circuit.ccx(c3, ce1, ce2)
    circuit.ccx(c1, c2, ce1)

def multi_controlled_x(control_qubits, aux_qubit, target_qubit, circuit):
    num_controls = len(control_qubits)
    
    for i in range(num_controls - 1):
        circuit.ccx(control_qubits[i], control_qubits[i + 1], aux_qubit)
    
    circuit.ccx(aux_qubit, control_qubits[-1], target_qubit)
    
    for i in range(num_controls - 2, -1, -1):
        circuit.ccx(control_qubits[i], control_qubits[i + 1], aux_qubit)
    


def Diff_with_auxiliary(qubits, auxiliary_qubits, circuit):
    for i in range(len(auxiliary_qubits)):
        circuit.ccx(qubits[i], qubits[i + 1], auxiliary_qubits[i])
    circuit.cz(qubits[-1], auxiliary_qubits[-1])
    for i in range(len(auxiliary_qubits) - 1, -1, -1):
        circuit.ccx(qubits[i], qubits[i + 1], auxiliary_qubits[i])

def Diff(*qubits, QuantumCircuit):
    circuit = QuantumCircuit
    for q in qubits:
        circuit.h(q)
        circuit.x(q)
    if len(qubits) == 4:
        CCCZ(qubits[0], qubits[1], qubits[2], qubits[3], circuit)
    else:
        auxiliary_qubits = [circuit.num_qubits - i - 1 for i in range(len(qubits) - 4)]
        Diff_with_auxiliary(qubits, auxiliary_qubits, circuit)
    for q in qubits:
        circuit.x(q)
        circuit.h(q)

def SizecountInv(QuantumCircuit):
    circuit = QuantumCircuit
    circuit.ccx(cl22,cr22,ce4)
    circuit.ccx(cl22,cr21,ce2)
    circuit.ccx(cl21,cr22,ce2)
    circuit.ccx(cl21,cr21,ce1)
    
    circuit.ccx(u2,cr11,cr22)
    circuit.ccx(u2,cr10,cr21)
    C_CX(u2,cr11,cr21,circuit)
    C_CX(u2,cr10,cr20,circuit)
    circuit.ccx(u1,bic,cr11)
    C_CX(u1,bic,cr10,circuit)
    
    circuit.ccx(v2,cl11,cl22)
    circuit.ccx(v2,cl10,cl21)
    C_CX(v2,cl11,cl21,circuit)
    C_CX(v2,cl10,cl20,circuit)
    circuit.ccx(v1,bic,cl11)
    C_CX(v1,bic,cl10,circuit)
    circuit = QuantumCircuit
               
        
def Bicheck(QuantumCircuit):
    circuit = QuantumCircuit

    circuit.ccx(v1,u1,e1)
    circuit.ccx(v2,u1,e3)
    circuit.ccx(v2,u2,e4)

    index = 1
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            circuit.ccx(globals()[f"v{i}"], globals()[f"u{j}"], globals()[f"e{index}_"])
            index += 1

    for i in range(1, n * m + 1):
        circuit.cx(globals()[f"e{i}"], globals()[f"e{i}_"])


    for i in virtues:
        circuit.x(i)    
    
    circuit.ccx(e1_,e2_,ce1)
    circuit.ccx(e3_,e4_,ce2)
    circuit.ccx(ce1,ce2,bic)
    circuit.ccx(e3_,e4_,ce2)
    circuit.ccx(e1_,e2_,ce1)
    
#     multi_controlled_x(virtues, ce1, bic, circuit)
    
    for i in virtues:
        circuit.x(i)



        
def Sizecount(QuantumCircuit):
    circuit = QuantumCircuit
    C_CX(v1,bic,cl10,circuit)
    circuit.ccx(v1,bic,cl11)
    C_CX(v2,cl10,cl20,circuit)
    C_CX(v2,cl11,cl21,circuit)
    circuit.ccx(v2,cl10,cl21)
    circuit.ccx(v2,cl11,cl22)

    C_CX(u1,bic,cr10,circuit)
    circuit.ccx(u1,bic,cr11)
    C_CX(u2,cr10,cr20,circuit)
    C_CX(u2,cr11,cr21,circuit)
    circuit.ccx(u2,cr10,cr21)
    circuit.ccx(u2,cr11,cr22)

    circuit.ccx(cl21,cr21,ce1)
    circuit.ccx(cl21,cr22,ce2)
    circuit.ccx(cl22,cr21,ce2)
    circuit.ccx(cl22,cr22,ce4)
    
    
def DiffInv(*qubits, QuantumCircuit):
    circuit = QuantumCircuit
    if len(qubits) < 4:
        raise ValueError("DiffInv requires at least 4 qubits as input.")
    reversed_qubits = list(reversed(qubits))
    for q in reversed_qubits:
        circuit.h(q)
        circuit.x(q)
    if len(reversed_qubits) == 4:
        if not callable(globals().get("CCCZ")):
            raise NameError("Function CCCZ is not defined.")
        globals()["CCCZ"](reversed_qubits[0], reversed_qubits[1], reversed_qubits[2], reversed_qubits[3], circuit)
    else:
        if not callable(globals().get("Diff_with_auxiliary")):
            raise NameError("Function Diff_with_auxiliary is not defined.")
        auxiliary_qubits = [circuit.num_qubits - i - 1 for i in range(len(reversed_qubits) - 4)]
        globals()["Diff_with_auxiliary"](reversed_qubits, auxiliary_qubits, circuit)
    for q in reversed_qubits:
        circuit.x(q)
        circuit.h(q)


def BicheckInv(QuantumCircuit):
    circuit = QuantumCircuit

    for i in virtues:
        circuit.x(i)

    circuit.ccx(e1_, e2_, ce1)
    circuit.ccx(e3_, e4_, ce2)
    circuit.ccx(ce1, ce2, bic)
    circuit.ccx(e3_, e4_, ce2)
    circuit.ccx(e1_, e2_, ce1)

    for i in virtues:
        circuit.x(i)


    for i in range(n * m, 0, -1):
        circuit.cx(globals()[f"e{i}"], globals()[f"e{i}_"])


    index = n * m
    for i in range(n, 0, -1):
        for j in range(m, 0, -1):
            circuit.ccx(globals()[f"v{i}"], globals()[f"u{j}"], globals()[f"e{index}_"])
            index -= 1


    circuit.ccx(v2, u2, e4)
    circuit.ccx(v2, u1, e3)
    circuit.ccx(v1, u1, e1)


    


def Oracle(size, QuantumCircuit):
    circuit = QuantumCircuit
#     Bicons(circuit)
    Bicheck(circuit)
    Sizecount(circuit)
    circuit.cx(size,o)
    SizecountInv(circuit)
#     BiconsInv(circuit)
    BicheckInv(circuit)
    


def Grover(size, *qubits, QuantumCircuit):
    circuit = QuantumCircuit
    if len(qubits) < 4:
        raise ValueError("Grover requires at least 4 qubits for Diffusion operation.")
    if not callable(globals().get("Oracle")):
        raise NameError("Oracle function is not defined in the current scope.")
    globals()["Oracle"](size, circuit)
    if not callable(globals().get("Diff")):
        raise NameError("Diff function is not defined in the current scope.")
    globals()["Diff"](*qubits[:4], QuantumCircuit=circuit)

    


# ## Circuit Construct

# In[2]:


################ simulator qubit count ################
# # Comment out this section during algorithm run
# qubit_num = 1
# max_qubits = 1000
# while qubit_num <= max_qubits:
#     try:
#         simulator = Aer.get_backend('qasm_simulator')
#         circuit = QuantumCircuit(qubit_num)
#         circuit.h(range(qubit_num))
#         circuit.measure_all()
#         transpiled_circuit = transpile(circuit, simulator)
#         result = simulator.run(transpiled_circuit, method='matrix_product_state').result()
#         if qubit_num % 10 == 0:
#             print(f"Simulated successfully with {qubit_num} qubits.")
#         qubit_num += 1
#     except Exception as e:
#         print(f"Simulation failed with {qubit_num} qubits. Error: {e}")
#         break
################ simulator qubit count ################

circuit = QuantumCircuit(1000,4)
        
# Bicheck(circuit)
# Sizecount(circuit)

# SizecountInv(circuit)
# BicheckInv(circuit)
# Oracle(ce1,circuit)
# Diff(v1,v2,u1,u2,circuit)

# Statepre(v1,v2,u1,u2,o,QuantumCircuit = circuit)
Statepre(*range(n*m),o,QuantumCircuit = circuit)

# Grover(ce1,*range(n*m),circuit)

# Grover(ce2,*range(n*m),QuantumCircuit = circuit)
# Grover(ce2,*range(n*m),QuantumCircuit = circuit)
for i_ in range(2):
    Grover(ce2, *range(n * m), QuantumCircuit=circuit)
    
# for i_ in range(3):
#     Grover(ce4,*range(n*m),circuit)

# Grover(ce4,*range(n*m),circuit)
# Grover(ce4,*range(n*m),circuit)
# Grover(ce4,*range(n*m),circuit)

circuit.measure([0,1,2,3],[0,1,2,3])


# In[3]:


# circuit.draw()


# In[4]:


# simulator = QasmSimulator()
# complied_circuit = transpile(circuit,simulator)
# job = simulator.run(complied_circuit, shots=1000)
# result = job.result()
# counts = result.get_counts(circuit)
# print(counts)

#############old version##############
# # IBMQ.save_account(TOKEN)
# IBMQ.load_account()
# IBMQ.providers()
# provider = IBMQ.get_provider(hub='ibm-q')
# # provider.backends()

# # backend = provider.get_backend('simulator_stabilizer')
# #error

# # backend = provider.get_backend('simulator_extended_stabilizer')
# # error

# backend = provider.get_backend('ibmq_qasm_simulator')

# # backend = provider.get_backend('simulator_mps')

# # backend = provider.get_backend('simulator_statevector')

# # backend = provider.get_backend('ibm_canberra')

# # backend.configuration()

# # qr = QuantumRegister(3)
# # cr = ClassicalRegister(3)
# # circuit = QuantumCircuit(qr, cr)
# # circuit.x(qr[0])
# # circuit.x(qr[1])
# # circuit.ccx(qr[0], qr[1], qr[2])
# # circuit.cx(qr[0], qr[1])
# # circuit.measure(qr, cr)
# # start = time.time()
# job = execute(circuit, backend, shots=20000)
#############old version##############


############## new version #############

# backend = Aer.get_backend('statevector_simulator')
# job = backend.run(circuit)

backend = Aer.get_backend('qasm_simulator')
simulation_options = {
    "method": "matrix_product_state"
}
job = backend.run(transpile(circuit, backend), shots=20000, **simulation_options)

############## new version #############



# In[5]:


# job.status()

result = job.result()
counts = result.get_counts()


# print(counts)

# end = time.time()
# print(end - start)

plot_histogram(counts)

# job.error_message()





