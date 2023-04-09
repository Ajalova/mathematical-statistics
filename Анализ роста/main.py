import numpy as np
import pandas as pd
import math
import openpyxl
import seaborn as sns
from scipy import special
from matplotlib import pyplot as plt
from scipy.stats import lognorm as normas
from scipy.stats import chi2

top_players = pd.read_excel('./asasa.xlsx')
top_players.head()
n=top_players.to_numpy()
a=np.array(n,float)
a=np.concatenate(a)

a_max=a.max()
a_min=a.min()
n=len(a)
print(n)
cont_of_intervals=int(1+math.log2(n))
long_og_intervals=(a_max-a_min)/cont_of_intervals
vibor_sred=np.sum(a)/n

#disp=np.var(a, ddof = 1)

s2=0
for i in range(0,n):
    s2=s2+(a[i]**2-vibor_sred)

s2=s2/(n+1)

b=[0]*(cont_of_intervals)
for i in range(0,n):
    e=0
    while(a[i]-a_min>(e+1)*long_og_intervals):
        e+=1
    b[e]+=1

mu=math.log(vibor_sred**2/math.sqrt(s2))
sigma=math.sqrt(math.log(s2/vibor_sred**2))

fig, ax = plt.subplots(figsize=(6, 5))
sns.histplot(a, stat='density', bins=cont_of_intervals)
x = [i for i in range(1, 120000)]
#plt.figure(2)
f1 = [1/(math.sqrt(2*math.pi)*sigma*(i))*(math.exp(1)**(-((math.log(i)-mu)**2)/(2*sigma**2))) for i in range(1, 120000)]

#for i in range(10000): f1.insert(0,0)
plt.plot(x, f1, ':b', label='1st component')

granici=[0]*(cont_of_intervals+1)
for i in range(0,cont_of_intervals+1):
    granici[i]=a_min+(i)*long_og_intervals

#plt.plot(granici, b, 'black')
def cdf(x):
    return 1/2+special.erf((math.log(x)-mu)/2/math.sqrt(sigma))/2

granici[cont_of_intervals]=(1000000)
granici[0]=-1000000

b0=[granici[i] for i in range(len(granici)-1)]
b1=[granici[i+1] for i in range(len(granici)-1)]
teor_thast=[n*(normas.cdf((b1[i])/math.exp(mu),sigma,0)-normas.cdf((b0[i])/math.exp(mu),sigma,0)) for i in range(1,len(granici)-1)]
teor_thast.insert(0,n*normas.cdf((b1[0])/math.exp(mu),sigma,0))


df = pd.DataFrame({'интервал': b0,'': b1,
                   'Имперические частоты': b,
                   'Теоретические частоты': teor_thast})
print(df)
print([((b[i]-teor_thast[i])**(2))/teor_thast[i] for i in range(len(teor_thast))])
ksi_kvadrat_vibor=sum([((b[i]-teor_thast[i])**(2))/teor_thast[i] for i in range(len(teor_thast))])
m=(cont_of_intervals-1)-2-1
print("chi^2_в=",ksi_kvadrat_vibor)
print("сhi^2(m)=",chi2.ppf(0.98,m))

plt.show()


