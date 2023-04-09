import numpy as np
import pandas as pd
import math
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import norm as normas
from scipy.stats import chi2

top_players = pd.read_excel('./asasa.xlsx')
top_players.head()
n=top_players.to_numpy()
a=np.array(n,float)
a=np.concatenate(a)
a_max=a.max()
a_min=a.min()
cont_of_intervals=int(1+math.log2(100))
long_og_intervals=(a_max-a_min)/cont_of_intervals
vibor_sred=np.sum(a)/100
#disp=np.var(a, ddof = 1)

s2=0
for i in range(0,100):
    s2=s2+(a[i]**2-vibor_sred)

s2=s2/101

b=[0]*8
for i in range(0,100):
    e=0
    while(a[i]-a_min>(e+1)*long_og_intervals):
        e+=1
    b[e]+=1

print(b)

#plt.figure(1)
fig, ax = plt.subplots(figsize=(6, 5))
sns.histplot(a, stat='density', bins=cont_of_intervals)
x = np.linspace(-10, 5)
#plt.figure(2)
f1 = 1/(math.sqrt(2*math.pi*s2))*(math.exp(1)**(-((x-vibor_sred)**2)/(2*s2)))
plt.plot(x, f1, ':b', label='1st component')


print(long_og_intervals)

#plt.figure(1)
granici=[0]*8
for i in range(0,8):
    granici[i]=a_min+(i)*long_og_intervals


#plt.plot(granici, b, 'black')

n=len(a)
granici.append(1000)
b.insert(0,0)
granici.insert(0,-1000)
b0=[granici[i] for i in range(len(granici)-1)]
b1=[granici[i+1] for i in range(len(granici)-1)]
teor_thast=[n*(normas.cdf((b1[i]-vibor_sred)/math.sqrt(s2))-normas.cdf((b0[i]-vibor_sred)/math.sqrt(s2))) for i in range(1,8)]
teor_thast.insert(0,n*normas.cdf(((b1[0]-vibor_sred)/math.sqrt(s2))))
teor_thast.append(n*(1-normas.cdf(((b0[8]-vibor_sred)/math.sqrt(s2)))))

df = pd.DataFrame({'интервал': b0,'': b1,
                   'Имперические частоты': b,
                   'Теоретические частоты': teor_thast})
print(df)

ksi_kvadrat_vibor=sum([(b[i]-teor_thast[i])**(2)/teor_thast[i] for i in range(len(teor_thast))])
m=(cont_of_intervals-1)-2-1
print("chi^2_в=",ksi_kvadrat_vibor)
print("сhi^2(m)=",chi2.ppf(0.98,m))

plt.show()


