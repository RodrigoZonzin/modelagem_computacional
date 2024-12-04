import numpy as np
import matplotlib.pyplot as plt

def plot(passo, x, y):
    #colors = np.random.rand(xsize)
    #plt.scatter(x, y, c=, alpha=0.5)
    fig = plt.figure(figsize=(9,7))
    plt.plot(x,y,linewidth=2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title("Random walk")
    plt.savefig("figs/y" + str(passo) + ".png", bbox_inches="tight")
    plt.close()
    #plt.show()

def Randomwalk1D(nsteps, xsize, dx, y, traces): 
    xvalues = np.arange(0,xsize,dx)

    plot(0, xvalues, y)

    for x in xvalues:
        traces[x][0] = y[x]
    
    for i in range (1,nsteps): # para cada passo de tempo 
        for x in xvalues: # para cada posicao x do espaco 
            for m in range(0,int(y[x])): # para cada molecula na posicao x 

                p = np.random.uniform(0,1)
                if p < 0.5: 
                    if (x < xsize - 1):
                        y[x+1] += 1
                        y[x] -= 1
                        traces[x][i] = y[x]
                        traces[x+1][i] = y[x+1]
                elif p < 1:
                    if (x > 0):
                        y[x-1] += 1
                        y[x] -= 1
                        traces[x-1][i] = y[x-1]
                        traces[x][i] = y[x]

        if (i % 10 == 0):
            plot(i, xvalues, y)    

    plt.figure(figsize=(9,7))
    times = range (0,nsteps)
    plt.plot(times,traces[45][:], color='r', label='y[45]')
    plt.plot(times,traces[70][:], color='m', label='y[70]')
    plt.xlabel('t')
    plt.ylabel('y')
    plt.legend(loc='best')
    plt.title("Evolution of y concentration at different points")
    plt.savefig("figs/y_x.png", format="png", bbox_inches="tight")


nsteps = 600 
xsize = 100 # tamanho do dominio (de 0 a 100)
dx = 1 # discretizacao Espacial 
npoints = int(xsize/dx)
num = 50 # numero de moleculas 

y = np.zeros(npoints)  # concentracao em cada ponto x do espaco   
y[int(xsize/2)] = num

traces = np.zeros((npoints,nsteps))

plt.figure(figsize=(9,7))
Randomwalk1D(nsteps,xsize,dx,y, traces) 
