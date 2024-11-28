import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import math
import random

pdf = PdfPages('sir-gillespiev1.pdf')
tf = 10
r = 0.3     #taxa de reproducao da presa
a = 0.05    #taxa de predacao
b = 0.01     #taxa de reproducao dos predadores
m = 0.09     #taxa de mortalidade dos predadores 

fig1 = plt.figure(figsize=(6,4))
plt.title('Predador')
plt.xlabel('tempo')
plt.ylabel('população')
ax1 = fig1.gca()

fig2 = plt.figure(figsize=(6,4))
plt.title('Presa')
plt.xlabel('tempo')
plt.ylabel('população')
ax2 = fig2.gca()

"""
fig3 = plt.figure(figsize=(6,4))
plt.title('Recuperados')
plt.xlabel('tempo')
plt.ylabel('população')
ax3 = fig3.gca()
"""

k = 0
colorindex = 0.
while k < 10:  
    random.seed(k)  
    times = []
    hresult = []
    presult = []

    h = 190
    p = 2
    i = 10
    t = 0
    hresult.append(h)
    presult.append(p)

    times.append(t)    
    while t < tf:

        R1 = r*h            #aumento da presa
        R2 =  a*h*p         #predacao
        R3 = b*h*p          #aumento do predador
        R4 = m*p            #morte do predador

        R = R1 + R2 +R3 + R4
        ran2 = random.uniform(0,1)

        #hunter
        if ran2 < R1/R:
            h += 1 

        elif ran2 < (R1+R2)/R:
            h -= 1

        elif ran2 < (R1+R2+R3)/R:
            p += 1

        else:
            p -+ 1


        hresult.append(h)
        presult.append(p)

        ran = random.uniform(0,1)
        tau = - math.log(ran)/R
        t = t + tau
        times.append(t)
        #print('t = ' + str(t) + '\n')

    ax1.plot(times,presult, color=(colorindex,colorindex,colorindex), label="S");
    ax2.plot(times,hresult, color=(colorindex,colorindex,colorindex), label="I")
    #ax3.plot(times,rresult, color=(colorindex,colorindex,colorindex), label="R")
    colorindex += 0.1
    k += 1

pdf.savefig(fig1)
pdf.savefig(fig2)
#pdf.savefig(fig3)
pdf.close()