from cProfile import label
import matplotlib.pyplot as plt
import os


with open('./results/michigan_nosel/exp3/convergence.txt') as f:
    lines = f.readlines()

with open('./results/michigan_10worst/exp1/convergence.txt') as f:
    lines2 = f.readlines()

y_GRIEF_DE = [int(x) for x in lines]

y_GRIEF_DE_sel = [int(x) for x in lines2]

list_files = os.listdir("/home/adriel/repos/GRIEF/grief_history1_michigan2")

list_files = [x.replace(".txt","") for x in list_files]

list_numeration = []
list_fit = []

count = 0
for x in list_files:
    list_numeration.append(x.split("_")[0])
    list_fit.append(x.split("_")[1])
    count+=1

aux = [int(x) for x in list_fit]
y_order = [int(x) for x in list_numeration]

y_GRIEF = [x for _, x in sorted(zip(y_order, aux))]

k = 9
plt.plot(y_GRIEF[k:1500+k], label='GRIEF')
plt.plot(y_GRIEF_DE, label = 'GRIEF-DE sem sel')
plt.plot(y_GRIEF_DE_sel, label = 'GRIEF-DE')
plt.legend()
plt.show()