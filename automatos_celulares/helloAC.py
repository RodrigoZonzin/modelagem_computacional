import numpy as np 
import matplotlib.pyplot as plt 

n = 3
t_size = 10
x = 0
y = 0

#celula morta 0 
#celula viva 1 

Ac = np.array([[1, 0, 0, 1], 
      [0, 1, 1, 1],
      [1, 1, 1, 0],
      [1, 1, 1, 1]])

Ac_next = np.copy(Ac)

def num_vizinhos_vivos(i, j): 
    k = 0 
    global Ac, n
    
    #vizinho da esquerda
    if Ac[i-1 % n][j] == 1: 
        k += 1
    
    #vizinho da direita
    if Ac[i+1  % n][j] == 1: 
        k += 1

    #vizinho de baixo
    if Ac[i][j-1  % n] == 1: 
        k += 1

    #vizinho de baixo
    if Ac[i][j+1  % n] == 1: 
        k += 1
    
    return k


for t in range(t_size): 
    for x in range(n):
        for y in range(n): 
            if Ac[x][y]  == 1 and num_vizinhos_vivos(x, y) < 2:
                Ac_next[x][y] = 0
            
            if Ac[x][y]  == 1 and num_vizinhos_vivos(x, y) > 3:
                Ac_next[x][y] = 0

            if Ac[x][y]  == 0 and num_vizinhos_vivos(x, y) == 3:
                Ac_next[x][y] = 1

            if Ac[x][y]  == 1 and num_vizinhos_vivos(x, y) == 3 or num_vizinhos_vivos(x, y) == 2:
                Ac_next[x][y] = 1

    print(Ac)
    print("++++++++"*5+'\n')
    plt.imshow(Ac, cmap = 'viridis')
    plt.savefig(f'ac_t_{t}.png')
    Ac= np.copy(Ac_next)




