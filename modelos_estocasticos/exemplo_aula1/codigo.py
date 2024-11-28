import numpy as np
import random 
import matplotlib.pyplot as plt



for k in range(5):
    y = 30
    yresult = []

    for t in range(50): 
        p = random.random()

        if p <= 0.5:
            y += 1
        else: 
            y -= 1
        
        yresult.append(y)


    yresult = np.array(yresult)

plt.scatter(y = yresult, x = range(50))
plt.ylim([20, 50])
plt.savefig(f'modelos_estocasticos/teste.png')
plt.plot(yresult)