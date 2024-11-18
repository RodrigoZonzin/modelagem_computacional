import numpy as np
import random 
import matplotlib.pyplot as plt


r = 0.5
m = 1-r

for k in range(5):
    y = 30
    tresult = []
    yresult = []

    for t in range(50): 
        p = random.random()
        if(y == 0): break
        p1 = r*y
        p2 = m*y

        s = p1+p2
        print(p, p1, p2, s)

        if p <= p1/s:
            y += 1
        else: 
            y -= 1
        
        yresult.append(y)
        tresult.append(t)


    yresult = np.array(yresult)

    plt.ylim([0, 50])
    plt.plot(yresult)

plt.savefig(f'modelos_estocasticos/testeCodigo2.png')
plt.show()
