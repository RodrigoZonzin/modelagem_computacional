import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import math
import random
import sys 

execucao_i = str(sys.argv[1])

pdf = PdfPages(f'sir-gillespiev1_teste{execucao_i}.pdf')
tf = 10
p = 0.1             #nascimentos
beta = 0.6          #transmissao 
alpha = 0.01        #taxa de resuscetibilidade  
gama = 0.01          #taxa de recuperacao
d = 0.5             #taxa de mortalidade devido à doenca

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
    i = 20
    r = 0

    t = 0

    sresult.append(s)
    iresult.append(i)
    rresult.append(r)
    times.append(t)    

    while t < tf:        
        N = s + i + r

        #nascimentos naturais
        R1 = p*N
        #infeccao 
        R2 = -beta*s*(i/N)
        #voltam a ser sucetiveis
        R3 = alpha*r      

        #infeccao
        R4 = beta*s*(i/N)
        #recuperacao
        R5 = -gama*i
        #morte ocasionada pela doenca
        R6 = -d*i

        #taxa de recuperacao
        R7 = gama*i 
        #voltam a ser sucetiveis
        R8 = -alpha*r

        R = R1 + R2 + R3 + R4 + R5 + R6 + R7 + R8

        ran2 = random.uniform(0,1)

        if ran2 < R1/R:
            s += 1
            
        elif ran2 < (R1+R2)/R:
            s -+ 1
            i += 1
            

        elif ran2 < (R1+R2+R3)/R:
            s += 1
            r -+ 1

        elif ran2 < (R1+R2+R3+R4)/R:
            i += 1
            s -= 1

        elif ran2 < (R1+R2+R3+R4+R5)/R:
            i -= 1
            r += 1
        
        elif ran2 < (R1+R2+R3+R4+R5+R6)/R:
            i -= 1

        elif ran2 < (R1+R2+R3+R4+R5+R6+R7)/R:
            r += 1
            i -= 1

        elif ran2 < (R1+R2+R3+R4+R5+R6+R7+R8)/R:
            r -+ 1
            s -+ 1

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