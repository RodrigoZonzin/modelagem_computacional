import numpy as np
import random 
import matplotlib.pyplot as plt

H = 20
P = 2 
r = 0.1
a = 0.1 
m = 0.2

for k in range(5):
    y = 30
    tresult = []
    hresult = []
    presult = []


    for t in range(50): 
        p = random.random()
        if(y == 0): break
        p1 = r*H
        p2 = a*H*p
        p3 = m*p
        s = p1+p2+p3

        
        if p <= p1/s:
            y += 1
        
        elif p<= (p1+p2)/s:
            H -= 1
            P += P

        else: 
            P -= 1
        
        hresult.append(H)
        presult.append(P)
        tresult.append(t)


    hresult = np.array(hresult)
    presult = np.array(presult)
    tresult = np.array(tresult)

    #plt.ylim([0, 50])
    plt.subplot(2, 2, 1)
    plt.plot(presult)
    plt.subplot(2, 2, 2)
    plt.plot(hresult)

plt.savefig(f'modelos_estocasticos/testePredadorPresa.png')
#plt.show()
