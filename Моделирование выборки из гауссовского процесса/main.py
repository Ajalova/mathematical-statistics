#!/usr/bin/env python
# coding: utf-8

# In[195]:


import pandas as pd
import numpy as np
from matplotlib import colors
import scipy.special as sc
import random
import math
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator


# данные
x = np.linspace(0, 5)
p = x*math.exp(1)**(-0.5*x*x)
n = 120



# ф-я распределения
f_r=1-math.exp(1)**(-0.5*x**2)

# обратная ф-я
def X(z):
    return math.sqrt(-math.log(1-z,math.e)*2)


# Моделируем вектор из n случайных чисел
y = [random.random() for i in range(n)]
print("вектор из n случайных чисел - ",[round(v,4) for v in y])

y1 = [X(y[i]) for i in range(n)]
print("преобразованный вектор - ",[round(v,4) for v in y1])

a=np.array(y1,float)


# Крайние члены вариационного ряда
a_max=a.max()
a_min=a.min()
print("max = ",a_max)
print("min = ",a_min)


# Число интервалов группировки по правилу Стёрджеса
cont_of_intervals=int(1+math.log2(n))
print("Число интервалов = ",cont_of_intervals)


# Длинна интервалов и выборочное среднее
long_intervals=(a_max-a_min)/cont_of_intervals
vibor_sred=np.sum(a)/n
print("Длинна интервалов = ",long_intervals)
print("выборочное среднее = ",vibor_sred)


# Дисперсия
disp=0
for i in range(0,n):
    disp=disp+(a[i]-vibor_sred)**2

disp=disp/(n-1)
print("Дисперсия = ",disp)


# частоты
b=[0]*cont_of_intervals
for i in range(0,n):
    e=0
    while(a[i]-a_min>(e+1)*long_intervals):
        e+=1
    b[e]+=1
print("Частоты- ",[v for v in b])
print(sum(b))


# гистограмма
plt.figure(1)
fig, ax = plt.subplots(figsize=(6, 5))
N, bins, patches=ax.hist(a, bins=cont_of_intervals)
fracs = ((N ** (1 / 5)) / N.max())
norm = colors.Normalize(fracs.min(), fracs.max())
for thisfrac, thispatch in zip(fracs, patches):
    color = plt.cm.viridis(norm(thisfrac))
    thispatch.set_facecolor(color)

# функция распределения
plt.figure(2)
plt.plot(x, p*100, ':b', label='1st component')


# Накопленные частоты

acc_rel_fq = [0]
for i in range (cont_of_intervals):
    acc_rel_fq.append(acc_rel_fq[i] + b[i]/n)
acc_rel_fq.pop(0)
print("Накопленные частоты- ",[round(v,4) for v in acc_rel_fq])



plt.figure(1)
plt.plot(x, f_r, ':b')
x = np.linspace(-1, 6)
def ind(z):
    ans=0
    if (z>0):ans=1
    return ans

def raspred_emper(z):
    ans=0
    for i in range(n):
        ans+=ind(z-a[i])/n
    return ans

x_param_e = np.linspace(-1, 6, 500)
y_param_e = np.vectorize(raspred_emper, otypes=[float])
graph_e = plt.plot(x_param_e, y_param_e(x_param_e), linewidth = 2)
ax = plt.gca()



g=0.1
epsilon=math.sqrt(1/2/n*math.log(2/g,math.e))
def r(x):
    e=raspred_emper(x)
    if (e+epsilon<1): return e+epsilon
    else: return 1
def l(x):
    e = raspred_emper(x)
    if (e-epsilon>0): return e-epsilon
    else: return 0
y_param_e = np.vectorize(r, otypes=[float])
graph_e = plt.plot(x_param_e, y_param_e(x_param_e), linewidth = 2)
ax = plt.gca()
y_param_e = np.vectorize(l, otypes=[float])
graph_e = plt.plot(x_param_e, y_param_e(x_param_e), linewidth = 2)
ax = plt.gca()

plt.show()
