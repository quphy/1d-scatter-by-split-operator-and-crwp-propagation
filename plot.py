import pandas as pd
import matplotlib.pyplot as plt

dfA = pd.read_csv('probability_sincdvr_crwp', sep='\s+', header=None, names=['x', 'y'])
dfB = pd.read_csv('probability_sincdvr_so', sep='\s+', header=None, names=['x', 'y'])

plt.figure(figsize=(8, 5))
plt.plot(dfA['x'], dfA['y'], label='crwp')
plt.plot(dfB['x'], dfB['y'], label='so')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('crwp vs so')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

