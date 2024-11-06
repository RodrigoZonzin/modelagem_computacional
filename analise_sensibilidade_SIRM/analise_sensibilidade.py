import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import scipy
from tqdm import tqdm
from uqtools import *
import SALib
from SALib.sample import saltelli
from SALib.analyze import sobol

#TEMPO
dt = 0.1
tfinal = 50.0
nsteps = int((tfinal)/dt)
t = np.arange(0, tfinal+dt, dt)

beta = 0.002 
alpha = 0.25 
l = 0.2 
d = 0.1
r = 0.1

"""
class Parameter:
    def __init__(self, name, value):
        self.name = name 
"""

def odeSystem(u, t, beta, alpha, l, d, r):
    global model_args
    S, I, R = u
    N = S + I + R 

    dSdt = -beta*S*I + l*R + r*N
    dIdt = beta*S*I - alpha*I - d*I 
    dRdt = alpha*I - l*R 

    return np.array([dSdt, dIdt, dRdt])

def eval_model(params):
    global model_args
    
    beta = params[0]
    alpha = params[1]
    l = 0.0 #params[2]
    d = params[3]
    r = 0.0 #params[4]

    model_args = (beta, alpha, l, d, r)

    S = 500
    I = 5
    R = 0
    u = (S, I, R)
    y = odeint(odeSystem, u, t, args=model_args)

    return y[:,0], y[:,1], y[:,2]

if __name__ == "__main__":
    opt_save_evals = True

    opt_sobol_salib = True    

    label_param = ['beta', 'alpha', 'lambda', 'd', 'r']

    # Sensitivity analysis
    if(opt_sobol_salib):
        min_bound = 0.1  
        max_bound = 10.0 

        #Valor de referência para gerar as amostras
        model_params = (beta, alpha, l, d, r)
        vbounds = [] #Vetor de limites inferior e superior pra cada parâmetro
        for i in range(len(model_params)):
            vbounds.append([model_params[i]*min_bound, model_params[i]*max_bound])

        #Pra conseguir gerar a amostra, precisa dessa estrutura
        Nvars = 5 #numero de parametros
        # Define the model inputs
        problem = {
            'num_vars': Nvars,
            'names': label_param, 
            'bounds': vbounds
        }

        # Generate samples
        nsobol = 512 # 1024 #numeros de amostras

        param_values = saltelli.sample(problem, nsobol, calc_second_order=False)

        # Run model (example)
        k = 0
        nexec = np.shape(param_values)[0] #numero de linhas da matriz com as amostras
        lpts = range( nexec )
        #Indice sobol principal
        
        sm_S = []
        sm_I = []
        sm_R = []

        print("evaluating samples for SA: ")
        #tqdm exibe a barra de progresso
        #Uma amostra é um conjunto de valores de parâmetros
        for i in tqdm(lpts,bar_format='{l_bar}{bar:20}{r_bar}{bar:-20b}'):
        #for s in samples.T:
            s = np.array(param_values[k,:]) #acessa todas as colunas da linha k
            k = k+1
            
            S, I, R = eval_model(s) #resolve as EDOs para a amostra na linha k
            
            #ignora 1o passo de tempo porque é a condição inicial e não tem variancia, levando a um NaN
            sm_S.append(S[1:])
            sm_I.append(I[1:])
            sm_R.append(R[1:])

        sm_S = np.array(sm_S)
        sm_I = np.array(sm_I)
        sm_R = np.array(sm_R)

        #Indice S1 para cada passo de tempo
        S1_S = np.zeros((nsteps,Nvars))
        S1_I = np.zeros((nsteps,Nvars))
        S1_R = np.zeros((nsteps,Nvars))

        # Perform analysis
        for i in range(nsteps): #para cada valor no tempo
            
            S1_S[i,:] = sobol.analyze(problem, sm_S[:,i], calc_second_order=False, parallel=True, print_to_console=False)['S1']
            S1_I[i,:] = sobol.analyze(problem, sm_I[:,i], calc_second_order=False, parallel=True, print_to_console=False)['S1']
            S1_R[i,:] = sobol.analyze(problem, sm_R[:,i], calc_second_order=False, parallel=True, print_to_console=False)['S1']
            
            print('step = ' + str(i))

        print("salvando arquivos Sobol")
        if(opt_save_evals):
            np.savetxt('ouput_sobol_S.txt',S1_S)
            np.savetxt('ouput_sobol_I.txt',S1_I)
            np.savetxt('ouput_sobol_R.txt',S1_R)
        
        plt.figure()
        plt.title("S")
        plot_sensitivity_mc(plt, t[1:], S1_S.T, label_param)
        plt.legend(bbox_to_anchor=(0.3, 1.00), ncol=3, fancybox=True, prop={'size': 8})
        plt.tight_layout()
        plt.savefig('output_sens_S.pdf')

        plt.figure()
        plt.title("I")
        plot_sensitivity_mc(plt, t[1:], S1_I.T, label_param)
        plt.legend(bbox_to_anchor=(0.3, 1.00), ncol=3, fancybox=True, prop={'size': 8})
        plt.tight_layout()
        plt.savefig('output_sens_I.pdf')

        plt.figure()
        plt.title("R")
        plot_sensitivity_mc(plt, t[1:], S1_R.T, label_param)
        plt.legend(bbox_to_anchor=(0.3, 1.00), ncol=3, fancybox=True, prop={'size': 8})
        plt.tight_layout()
        plt.savefig('output_sens_R.pdf')

# Fim
print('done')