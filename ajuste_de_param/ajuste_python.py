import numpy as np
import pandas as pd
from scipy.integrate import  solve_ivp
from scipy.optimize import differential_evolution
import math

path = ''

def odeSystem(t, alpha, beta):
    beta, alpha, = constants
    S, I, R, = u

    N = S + I + R 

    dSdt = -beta*S*I 
    dIdt = beta*S*I - alpha*I 
    dRdt = alpha*I 

    return np.array([dSdt,dIdt,dRdt])
    

def isReferenceTime(times, ct):
    for t in times: 
        if (abs(ct - t) <= 10**(-5)):
            return True 
    return False

def solve(x):
    global data, reference_times
    dt = 0.01
    tfinal = 50
    times = np.arange(0,tfinal+dt,dt)

    N0 = x[0]
    u = [N0]

    r = x[1]
    k = x[2]
    params = (r, k)
    
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
            #sumobs = 1
            i = i + 1
        j = j + 1

    error = math.sqrt(error/sumobs) #Erro norma 2                  
    return error 

if __name__ == "__main__":
    global data, reference_times 
    data = pd.read_csv(path+'I.csv', delimiter=',')
    
    reference_times = data["Semana"]
    print(reference_times)
    bounds = [
        (0.001, 0.01),  #limites inferior e superior alpha
        (0.01, 1.0)     #limites inferior e superior alpha
    ]

    #chama evolução diferencial, result contém o melhor individuo
    solucao = differential_evolution(solve, bounds, strategy='rand2bin', maxiter=50, popsize=80,
        atol=10**(-3), tol=10**(-3), mutation=0.8, recombination=0.5, disp=True, workers=4)
    
    print(solucao.x)
    print(solucao.success)
    #saving the best offspring...
    np.savetxt('solucao_ajuste.txt',solucao.x, fmt='%.2f')        
    best=solucao.x
    error = solve(best)
    print("ERROR ", error)