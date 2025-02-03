try:
    from ca import *

except ModuleNotFoundError:
    from sys import path
    path.insert(0, '..')
    from ca import *
    from numpy.random import random


class ForestProp(CA):
    global QUEIMADO, QUEIMANDO, VAZIO, FLORESTA
    def rule(self, x, y):
        #estado atual
        s = self[x, y]

        #vizinhos
        n = neighbors8(self, x, y)
        
        #Estado queimado ou vazio: nao muda de estado.
        if s == VAZIO or s == QUEIMADO: 
            pass

        #Estado floresta: pode mudar para o estado queimando com probabilidade pf 
        #e se tiver pelo menos um vizinho queimando;
        if s == FLORESTA: 
            p = random()

            if p < pF and QUEIMANDO in n: 
                s == QUEIMANDO
        
        #Estado queimando: muda para o estado queimado ap´os um passo de tempo
        if s == QUEIMANDO: 
            s = QUEIMADO



VAZIO       = 0
FLORESTA    = 1
QUEIMANDO   = 2
QUEIMADO    = 3
pF          = 0.75
d           = 300                    #dimensao do dominio

nVazio      = int((1-pF)*d - 1)     #numero de celulas vazias = (1 - P[ser floresta])*d
nFloresta   = int(pF*d)


Ac = ForestProp(20, [VAZIO] + [FLORESTA]*0 + [QUEIMANDO]*0 + [QUEIMADO]*0)

i = 0 
while i < d: 
    
    Ac.addrandomvalues(VAZIO, FLORESTA, 300)
    i += 1


#plot(Ac, N=20, out="meuca.pdf")
plot(Ac, N=50, out='forestProp.pdf', graphic = True, vmax = 3,
    colors=['white', 'green', 'orange', 'black'],
    names=['VAZIO', 'FLORESTA', 'QUEIMANDO', 'QUEIMADO']
)