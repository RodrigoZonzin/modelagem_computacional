try:
    from ca import *
    from random import * 

except ModuleNotFoundError:
    from sys import path
    path.insert(0, '..')
    from ca import *
    from numpy.random import random


class ForestProp(CA):
    def rule(self, x, y):
        global QUEIMADO, QUEIMANDO, VAZIO, FLORESTA, pF

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

N           = 2
pF          = 0.45    
pFloresta   = 0.99
d           = int(pFloresta*N**2) 
                

Ac = ForestProp(N, values = 4, random_values=False)

#for i in range(d+200):
#    Ac.add(FLORESTA, points=[(randint(0, N), randint(0, N))], size = (1,1))

Ac.add(0, points = [(0, 0)], size = (1, 1))
Ac.add(1, points = [(0, 1)], size = (1, 1))
Ac.add(2, points = [(1, 0)], size = (1, 1))
Ac.add(0, points = [(1, 1)], size = (1, 1))
#xy_meio = int(N/2)
#print(xy_meio)
#Ac.add(QUEIMANDO, points=[(xy_meio, xy_meio)], size= (1,1))


plot(Ac, N=50, out='forestProp.pdf', graphic = True, vmax = 3, 
    colors = ['white', 'darkgreen', 'orange', 'black'],
    names = ['Vazio', 'Floresta', 'Queimando', 'Queimado'])
