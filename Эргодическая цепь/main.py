import numpy as np
import pandas as pd
import networkx as nx
from matplotlib import pyplot as plt
from numpy import linalg as LA
from scipy.linalg import solve
import random


#чтение таблицы в массив
m=4
k=20
N=150
top_players = pd.read_excel('./table2.xlsx')
top_players.head()
n=top_players.to_numpy()
a=np.array(np.delete(n, 0, 1),float)
print("матрица переходных вероятностей:")
print('\n'.join('\t'.join(map(str, row)) for row in a))

#печать графа

G1 = nx.MultiDiGraph()
G = nx.Graph()
for i in range(len(a)):
    for j in range(len(a)):
        if (a[i][j]!=0 and a[i][j]!=None):
            G.add_edge(i+1, j+1, weight=a[i][j])
            G1.add_edge(i + 1, j + 1, weight=a[i][j])

pos = nx.circular_layout(G1)
nx.draw_networkx_nodes(G, pos)
pos1={}
j=(0,0,{"":0})
pos1 = {}
pos2 = {}
for i in G.edges(None, True):
    if (i[0] == j[1] and j[0] == i[1]):
        pos1[i[0]] = pos[i[0]].copy()
        pos1[i[0]][0] += 0.255
        pos1[j[0]][0] -= 0.255
        pos2[i[0]] = (i[2]['weight'])
    if (i[0] == i[1]):
        pos1[i[0]] = pos[i[0]].copy()
        pos1[i[0]][1] += 0.255
        pos2[i[0]] = (i[2]['weight'])
    j=i


nx.draw_networkx_edges(G1, pos, edgelist=None, width=1,connectionstyle='arc3,rad=0.2')
edge_labels = nx.get_edge_attributes(G, "weight")
nx.draw_networkx_edge_labels(G, pos, edge_labels,font_size=8)

nx.draw_networkx_labels(G1, pos,font_size=10)
nx.draw_networkx_labels(G, pos1, {n:lab for n,lab in pos2.items() if n in pos1},font_size=8, font_color='k', font_family='sans-serif', font_weight='normal', alpha=None,  horizontalalignment='center', verticalalignment='bottom', ax=None, clip_on=False)
ax = plt.gca()
ax.margins(0.08)
plt.axis("off")
plt.tight_layout()

#эргодичность

P_2 = np.linalg.matrix_power(a, 2)
print("n = 2:\n", P_2)
P_3 = np.linalg.matrix_power(a, 3)
print("\nn = 3:\n", P_3)
P_500 = np.linalg.matrix_power(a, 20)
print("\nn = 20:\n", P_500)
P_1000 = np.linalg.matrix_power(a, 500)
print("\nn = 500:\n", P_1000)
P_1000 = np.linalg.matrix_power(a, 1000)
print("\nn = 1000:\n", P_1000)


# поиск предельных вероятностей

M=[]
for i in range(m):
    M.append([0]*(m))
    for j in range(m):
        if(i==j):
            M[i][j]=-1+a[j][i]
        else:
            M[i][j] = a[j][i]
M.append([1]*(m))
del M[m-1]
print("матрица для СЛАУ:")
print(np.matrix(M))
pi=solve(M, [0]*(m-1)+[1])
print("\nвектор предельных вероятностей:")
print('pi_j=',pi)

#моделирование начального вектора вероятностей состояний

r=sorted(np.random.uniform(low = 0.0, high = 1.0, size = m))
print('r=',r)
p0=[r[0]]
for i in range(1,m-1):
    p0.append(r[i]-r[i-1])
p0.append(1-r[m-2])
print("\nначальный вектор вероятностей:")
print('p_0=',p0)

#безусловные вероятности состояний смоделированной цепи на k шаге
Pk=LA.matrix_power (a.transpose(), k)
pk=np.matmul(Pk, p0)
print("\nбезусловные вероятности состояний смоделированной цепи на k шаге:")
print('p(',k,')=',pk)

#моделированиe траекторий
ss=[]
for i in range(N):
    s=[0]*k
    otr=np.cumsum(p0.copy())
    otr1=p0.copy()
    for j in range(k):
        r0 = random.uniform(0, otr[m-1])
        num=0
        if (r0 < otr[0]):
            num=0
        else:
            for i in range(1,m):
                if (otr[i-1]<r0<otr[i]):
                    num=i
        if(len(np.where(otr1==0.5)[0])==2):
            s[j]=1+random.choice(np.where(otr1==0.5)[0])
        else:
            s[j]=1+num
        otr = np.cumsum(a[s[j]-1])
        otr1=a[s[j]-1]

    print("S=",s)
    ss.append(s)
emp=[]
for i in range(k):
    emp.append([0]*m)
    for j in range(m):
        for l in range(N):
            if(ss[l][i]==j+1):
                emp[i][j]+=1/N

print('\n\n')
print('Эмперические вероятности')
print(np.matrix(emp))

fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot()
x = np.arange(m)
y3 = pi
y2 = emp[k-1]
y1 = pk
w = 0.3

fig.set_facecolor('floralwhite')
ax.bar(x +w, y1, width=w, label = "Теоретические вероятности")
ax.bar(x , y2, width=w, label = "Эмпирические вероятности")
ax.bar(x - w, y3, width=w, label = "Предельные вероятности")
ax.legend()
plt.show()

plt.show()





