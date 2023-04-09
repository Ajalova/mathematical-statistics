from scipy.stats import chi2
from scipy.stats import norm as normas
from scipy.special import stdtrit
import seaborn as sns
import numpy as np
import pandas as pd
import math
import openpyxl
from matplotlib import pyplot as plt

var1=[0.1,1.8,3,3,2.4,0.15]
primer=[0.05,-3.5,2,2,-4,0.1,-3.8]
t=primer
#уровень доверия
alpha=t[0]
#гипотеза Н0 а=а0 против H2:a>a0
a0=t[1]
#гипотеза Н01 σ=σ0 против H3:σ>σ0
sigma0=t[2]
#гипотеза Н02 а=а0 против H4:a<-4 при σ=2
sigma=t[3]
a1=t[4]
#
eps=t[5]

c1=t[6]



if (t==var1): top_players = pd.read_excel('./table.xlsx')
else: top_players = pd.read_excel('./table1.xlsx')

top_players.head()
n=top_players.to_numpy()
a=np.array(n,float)
a=np.concatenate(a)
count=len(a)

#Крайние члены вариационного ряда и размах выборки
a_max=a.max()
a_min=a.min()
cont_of_intervals=int(1+math.log2(count))
long_og_intervals=(a_max-a_min)/cont_of_intervals
# Выборочное среднее
vibor_sred=sum(a)/count
# Cреднее квадратичное отклонение
S2=sum([ (a[i]-vibor_sred)**2 for i in range(count)])/(count-1)
print('Выборочное среднее ',vibor_sred)
print('Cреднее квадратичное отклонение ',S2)

v_sigma=0
for i in range(0,count):
    V_sigma=v_sigma+(a[i]**2-vibor_sred)

v_sigma=v_sigma/count

v=[0]*cont_of_intervals
for i in range(0,count):
    e=0
    while(a[i]-a_min>(e+1)*long_og_intervals):
        e+=1
    v[e]+=1

#Относительные частоты
ot_v=[v[i]/count for i in range(cont_of_intervals)]
print('Сумма относительных частот',sum(ot_v))
#Плотность относительных частот
plotnost=[ f"{v[i]/count/long_og_intervals:.4f}" for i in range(cont_of_intervals)]

df = pd.DataFrame({'Интервал': [f"{a_min+long_og_intervals*(i):.3}" for i in range(0,cont_of_intervals)],
                   ' ': [f"{a_min+long_og_intervals*(i+1):.3}" for i in range(0,cont_of_intervals)],
                   'Частота': [v[i] for i in range(cont_of_intervals)],
                   'Относительная частота': [f"{ot_v[i]:.3}" for i in range(cont_of_intervals)],
                   'Плотность': plotnost}, index = [''] * (cont_of_intervals))
print(df)

plt.figure(1)
fig, ax = plt.subplots(figsize=(6, 5))
sns.histplot(a, stat='density', bins=cont_of_intervals)

#f1 = 1/(2*math.sqrt(math.pi*sigma))*(math.exp(1)**(-((x-vibor_sred)**2)/(2*sigma)))

C2=a0+stdtrit(count-1,alpha)/((count/S2)**(1/2))
print('Критическое множество - x̅>',C2)
if(vibor_sred>C2): print('Гипотеза Н0: а=',a0,' отклоняется, так как выборочное среднее принадлежит критическому множеству')
else: print('Гипотеза Н0: а=',a0,' принимается, так как выборочное среднее не принадлежит критическому множеству')

C3=chi2.ppf(alpha,count-1)*sigma0**2/(count-1)
print('Критическое множество - S^2>',C3)
if(S2>C2): print('Гипотеза Н0: σ=',sigma0,' отклоняется, так как cреднее квадратичное отклонение принадлежит критическому множеству')
else: print('Гипотеза Н0: σ=',sigma0,' принимается, так как реднее квадратичное отклонение не принадлежит критическому множеству')


C4=a0+normas.ppf(alpha)*sigma/(count)**(1/2)
print('Критическое множество - x̅>',C4)
if(vibor_sred>C2): print('Гипотеза Н0: а=',a0,' отклоняется, так как выборочное среднее принадлежит критическому множеству')
else: print('Гипотеза Н0: а=',a0,' принимается, так как выборочное среднее не принадлежит критическому множеству')

beta=1-normas.cdf((C4-a1)/sigma*math.sqrt(count))
print('β=',beta)

#Оптимальное значение a1', при котором ошибка второго рода не превышает ε
a1_opt=C4-normas.ppf(1-eps)/math.sqrt(count)*sigma
print('Оптимальное значение a1 при котором ошибка второго рода не превышает ε',a1_opt)

x = np.linspace(-10, 10)
f1 = 1/(2*math.sqrt(math.pi*sigma))*(math.exp(1)**(-((x-a0)**2)/(2*sigma)))
plt.plot(x, f1, ':b', label='N(a0,σ)')
f2 = 1/(2*math.sqrt(math.pi*sigma))*(math.exp(1)**(-((x-a1)**2)/(2*sigma)))
plt.plot(x, f2, ':r', label='N(a1,σ)')
plt.legend()

# Построим последовательный критерий Вальда для проверки гипотезы H_0: a = -3.5 = a_0
# против альтернативы H_1: a = -4 = a_1 при известном σ = 2 = σ_1.

A=(1-beta)/alpha
B=beta/(1-alpha)
def l_na_l(i):
     ans=((a0**(2)-a1**(2))/(2*sigma**(2))+(a1-a0)*a[0]/(sigma**(2)))
     for j in range(1,i+1):
         ans=ans+((a0**(2)-a1**(2))/(2*sigma**(2))+(a1-a0)*a[j]/(sigma**(2)))
     return math.exp(ans)

plt.figure(1)
x = [i for i in range(1,count+1)]
y= [l_na_l(i) for i in range(count)]
plt.plot(x, y, 'b', label='1st component')
plt.plot([0,120], [A,A],  label='1st component')
plt.plot([0,120], [B,B],  label='1st component')
c1=c1/count
print("математическое ожидание момента принятия решения при основной гипотезе H0", -(alpha*math.log(A)+(1-alpha)*math.log(B))/(a1-a0)**(2)*2*sigma**2)
print("математическое ожидание момента принятия решения при основной гипотезе H1", round((beta*math.log(B)+(1-beta)*math.log(A))/(a1-a0)**(2)*2*sigma**2,7))
c=math.exp(count*((a1-a0)*c1/sigma**(2)+(a0**(2)-a1**(2))/(2*sigma**(2)*count)))
print("критическое множество: S={L(X,a1)/L(X,a0)>=", round(c,5),"}")
plt.plot([0,120], [c,c],  label='1st component')

plt.show()
