import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import math
import random

pdf = PdfPages('sir-gillespiev1.pdf')
tf = 10
rs = 0.3 #nascimentos
ms = 0. #morte natural
beta = 0.01
nu = 0.2

fig1 = plt.figure(figsize=(6,4))
plt.title('Suscetíveis')
plt.xlabel('tempo')
plt.ylabel('população')
ax1 = fig1.gca()

fig2 = plt.figure(figsize=(6,4))
plt.title('Infectados')
plt.xlabel('tempo')
plt.ylabel('população')
ax2 = fig2.gca()

fig3 = plt.figure(figsize=(6,4))
plt.title('Recuperados')
plt.xlabel('tempo')
plt.ylabel('população')
ax3 = fig3.gca()


k = 0
colorindex = 0.
while k < 10:  
    random.seed(k)  
    times = []
    sresult = []
    iresult = []
    rresult = []
    s = 190
    i = 10
    r = 0
    t = 0
    sresult.append(s)
    iresult.append(i)
    rresult.append(r)
    times.append(t)    
    while t < tf:
        if i == 0 or s ==0: 
            break;        
        R1 = beta * s * i
        R2 = nu * i
        R3 = rs * s
        R4 = ms * s
        R = R1 + R2 + R3 + R4
        ran2 = random.uniform(0,1)

        if ran2 < R1/R:
            s = s - 1
            i = i + 1
        elif ran2 < (R1+R2)/R:
            i = i - 1
            r = r + 1
        elif ran2 < (R1+R2+R3)/R:
            s = s + 1
        elif ran2 < (R1+R2+R3+R4)/R:
            s = s - 1
        sresult.append(s)
        iresult.append(i)
        rresult.append(r)
        ran = random.uniform(0,1)
        tau = - math.log(ran)/R
        t = t + tau
        times.append(t)
        #print('t = ' + str(t) + '\n')

    ax1.plot(times,sresult, color=(colorindex,colorindex,colorindex), label="S");
    ax2.plot(times,iresult, color=(colorindex,colorindex,colorindex), label="I")
    ax3.plot(times,rresult, color=(colorindex,colorindex,colorindex), label="R")
    colorindex += 0.1
    k += 1

pdf.savefig(fig1)
pdf.savefig(fig2)
pdf.savefig(fig3)
pdf.close()