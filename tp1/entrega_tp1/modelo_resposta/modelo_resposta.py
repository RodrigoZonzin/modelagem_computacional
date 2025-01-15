# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def ode_system(t, u, constants):
    """
    system of first order differential equations
    t: discrete time step value
    y: state vector (vetor das variaveis/populacoes do modelo)
    """
    alpha   = constants[0] 
    mu_reg  = constants[1]
    beta    = constants[2]
    gama    = constants[3]
    mu_a    = constants[4]
    rho     = constants[5]
    alpha_a = constants[6]
    etha_ch = constants[7] 
    v       = constants[8] 
    etha_M  = constants[9] 
    w_reg   = constants[10] 
    etha_a  = constants[11] 

    TD  = u[0]
    N   = u[1]
    CH  = u[2]
    A   = u[3]
    M   = u[4]
    
    dTDdt   = alpha*N - mu_reg*M
    dNdt    = beta*TD + ((gama*CH)/(1+mu_a*A)) - alpha*N
    dCHdt   = (rho*N)/(1+mu_a*A) - etha_ch*CH
    dMdt    = v*N - etha_M *M
    dAdt    = w_reg*M - etha_a*A

    return np.array([dTDdt, dNdt, dCHdt, dMdt, dAdt])

def rk4(f, tk, _uk, _dt=0.01, **kwargs):
    """
    single-step fourth-order numerical integration (RK4) method
    f: system of first order ODEs
    tk: current time step
    _uk: current state vector [y1, y2, y3, ...]
    _dt: discrete time step size
    **kwargs: additional parameters for ODE system
    returns: y evaluated at time k+1
    """
    # evaluate derivative at several stages within time interval
    k1 = f(tk, _uk, **kwargs)
    k2 = f(tk + _dt / 2, _uk + (k1 * (_dt / 2)), **kwargs)
    k3 = f(tk + _dt / 2, _uk + (k2 * (_dt / 2)), **kwargs)
    k4 = f(tk + _dt, _uk + (k3 * _dt), **kwargs)

    # return an average of the derivative over tk, tk + dt
    return _uk + (_dt / 6) * (k1 + (2 * k2) + (2 * k3) + k4)

def euler(f, tk, _uk, _dt=0.001, **kwargs):
    """
    single-step explicit method
    func: system of first order ODEs
    tk: current time step
    _uk: current state vector [y1, y2, y3, ...]
    _dt: discrete time step size
    **kwargs: additional parameters for ODE system
    returns: y evaluated at time k+1
    """
    return _uk + _dt*f(tk, _uk, **kwargs)

def solveSystem(time, dt, y0, method):
    t = 0
    yk = [12, 3, 4, 1, 1]
    state = []
    """
    parameters = [0.95, 0.1 , 0.1 , 0.01, 
                  0.33, 0.92, 0.4, 0.8, 
                  0.1 , 0.3 , 0.01, 0.1 ]
    """

    global parameters
    print(parameters)

    if method == "euler":
        for t in time:
            state.append(yk)
            yk = euler(ode_system, t, yk, dt, constants=parameters)

    elif method == "rk4":
        for t in time:
            state.append(yk)
            yk = rk4(ode_system, t, yk, dt, constants=parameters)

    return np.array(state)

def save(time, state, names, filename):
    df = pd.DataFrame(state, columns = names)
    df.insert(0, 'time', time)
    df.to_csv(f'results{filename}.csv', float_format='%.5f', sep=',')

def plot(time, state, names, filename):
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 6)
    for i in range(len(names)):
        ax.plot(time, state[:, i], label=names[i], linewidth='2')

    #time2 = np.arange(0, tfinal + 0.001, 0.001)
    #sol = np.exp(r*time2)
    #ax.plot(time2, sol, label='solucao analitica', linewidth='2')
    ax.set(xlabel='Dias', ylabel='Populacao')
    plt.legend(loc='best')
    fig.savefig('modelo_inicial'+filename, format='png')
    plt.show()




if __name__ == "__main__":
    names = ['TD', 'N', 'CH', 'A', 'M']
    
    dt = 0.01
    tfinal = 50
    time = np.arange(0, tfinal + dt, dt)
    initial_condition = np.array([])

    parameters = np.random.rand((12))       #[0.1 0.2 0.4 0.5]
    nom_params = '_'.join(parameters)       #0.1_0.2_0.4_0.5

    result = solveSystem(time, dt, initial_condition, "rk4")
    

    save(time, result, names, nom_params)
    plot(time, result, names, nom_params)

