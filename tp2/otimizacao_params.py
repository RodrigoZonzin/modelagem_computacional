import numpy as np
import pandas as pd
from scipy.integrate import  solve_ivp
from scipy.optimize import differential_evolution
import math
import matplotlib.pyplot as plt 
import sys
import os

path = 'data/'

dt = 0.01
tfinal = 50
times = np.arange(0,tfinal+dt,dt)

S0 = 90082.0
I0 = 203.0 #203
R0 = 0.0

def odeSystem(t, u, beta, alpha):

    S, I, R = u
    dS_dt = - beta*S*I 
    dI_dt = beta*S*I - alpha*I
    dR_dt = alpha*I

    return [dS_dt, dI_dt, dR_dt]  

def isReferenceTime(times, ct):
    for t in times: 
        if (abs(ct - t) <= 10**(-5)):
            return True 
    return False

def solve(x):
    global data, reference_times

    u = [S0, I0, R0]

    beta = x[0]
    alpha = x[1]
    params = (beta, alpha)
    
    def solveOde(t, y):
        return odeSystem(t, y, *params)

    results = solve_ivp(solveOde,(0, tfinal), u, t_eval=times, method='RK45')    

    u = results.y[1,:]
    error = 0
    sumobs = 0
    i = 0
    j = 0    
    for t in times:
        if isReferenceTime(reference_times,t):
            p_data = data["I"][i]
            error += (u[j] - p_data)*(u[j] - p_data) 
            sumobs += p_data*p_data
            i = i + 1
        j = j + 1

    error = math.sqrt(error/sumobs) #Erro norma 2                  
    return error 

if __name__ == "__main__":
    global data, reference_times 
    
    beta_piso = 0.00001
    alpha_piso = 0.01

    data = pd.read_csv(path+'I.csv', delimiter=',')
    reference_times = data["Semana"]
    dados_I = data["I"]

    results = []

    for i in range(1, 4):
        for j in range (1, 4): 
            bounds = [(beta_piso*i, 0.01), (alpha_piso*j, 0.9)]



            #chama evolução diferencial, result contém o melhor individuo
            solucao = differential_evolution(solve, bounds, strategy='rand2bin', maxiter=50, popsize=40,atol=10**(-3), tol=10**(-3), mutation=0.8, recombination=0.5, disp=True, workers=6)


            best = solucao.x
            error = solve(best)

            results.append({'bounds': bounds,
                            'solucao': solucao.x,
                            'erro': error})

            print(f"{bounds}: {solucao.x}. Erro: {error}")

            u = [S0, I0, R0]
    
            df = pd.DataFrame(results)
            df.to_csv('resultados.csv')
        