import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys 
import os
import glob
import pandas as pd

data_path = 'results/'
fig_path = 'figs/'

x = np.loadtxt(data_path + 'x.csv')
#y = np.loadtxt(data_path + 'y.csv')
times = np.loadtxt(data_path + 't.csv')

values = [("N", "green"), ("B", "red"), ("Ch", "pink")]

for name, color in values: 
    file_name = name + "*.csv"
    all_files = glob.glob(os.path.join(data_path, file_name))
    
    for f in all_files:
        t = f.split("_")[1].split(".csv")[0]
        table = pd.read_csv(f)
        #print(table)
        fig_name = f.split("_")[0].split("/")[1]

        fig = plt.figure(figsize=(10,7))
        plt.title(fig_name)
        plt.plot(table["x"],table["value"], linewidth=2, color=color)
        plt.xlabel('x (mm)')
        plt.ylabel('y (concentration)')
        fig.savefig(fig_path + fig_name + '_' + t + '.png', format='png', bbox_inches='tight')
        plt.close()