import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('I.csv')

df['I'].plot()
plt.xlabel("Semana")
plt.ylabel("Populacao")
plt.savefig('dados.png')