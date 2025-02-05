try:
    from ca import *
    from random import * 

except ModuleNotFoundError:
    from sys import path
    path.insert(0, '..')
    from ca import *
    from numpy.random import random


class ForestProp(CA):
    global QUEIMADO, QUEIMANDO, VAZIO, FLORESTA, pF
    def rule(self, x, y):
        #estado atual
        s = self[x, y]

        #vizinhos
        n = neighbors8(self, x, y)
        
        #Estado queimado ou vazio: nao muda de estado.
        if s == VAZIO or s == QUEIMADO: 
            return s

        #Estado floresta: pode mudar para o estado queimando com probabilidade pf 
        #e se tiver pelo menos um vizinho queimando;
        if s == FLORESTA: 
            p = random()

            if p < pF and QUEIMANDO in n: 
                return QUEIMANDO
            
            return s
        
        #Estado queimando: muda para o estado queimado ap´os um passo de tempo
        if s == QUEIMANDO: 
            return QUEIMADO



VAZIO       = 0
FLORESTA    = 1
QUEIMANDO   = 2
QUEIMADO    = 3

N           = 50
pF          = 0.99   
d           = int(pF*N**2) 
                

Ac = ForestProp(N, values = 4, random_values=False)

for i in range(d):
    Ac.add(FLORESTA, points=[(randint(0, N), randint(0, N))], size = (1,1))

for i in range(1):
    Ac.add(VAZIO, points=[(randint(0, N), randint(0, N))], size = (1,1))

#for i in range(2): 
Ac.add(QUEIMANDO, points=[(0,0)], size = (1,1))
Ac.add(QUEIMANDO, points=[(0,1)], size = (1,1))

for i in range(5): 
    Ac.add(QUEIMADO, points=[(randint(0, N), randint(0, N))], size = (1,1))

plot(Ac, N=50, out='forestProp.pdf', graphic = True, vmax = 3, 
    colors = ['white', 'darkgreen', 'orange', 'black'],
    names = ['Vazio', 'Floresta', 'Queimando', 'Queimado'])