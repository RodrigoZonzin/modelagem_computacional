import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('casos_2021.csv')

df['I'].plot(kind = 'scatter')
plt.xlabel("Semana")
plt.ylabel("Populacao")
plt.savefig('dados_2021.png')