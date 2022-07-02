from cProfile import label
import matplotlib.pyplot as plt
import os
import re
from operator import add
from numpy import empty

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp1/convergence.txt') as f:
    lines = f.readlines()

y_GRIEF_DE_NOSEL1 = [int(x) for x in lines]

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp2/convergence.txt') as f:
    lines = f.readlines()

y_GRIEF_DE_NOSEL2 = [int(x) for x in lines]

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp3/convergence.txt') as f:
    lines = f.readlines()

y_GRIEF_DE_NOSEL3 = [int(x) for x in lines]

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp4/convergence.txt') as f:
    lines = f.readlines()

y_GRIEF_DE_NOSEL4 = [int(x) for x in lines]

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp5/convergence.txt') as f:
    lines = f.readlines()

y_GRIEF_DE_NOSEL5 = [int(x) for x in lines]

grief_data = list( map(add, y_GRIEF_DE_NOSEL1, y_GRIEF_DE_NOSEL2))
grief_data = list( map(add, grief_data, y_GRIEF_DE_NOSEL3))
grief_data = list( map(add, grief_data, y_GRIEF_DE_NOSEL4))
grief_data = list( map(add, grief_data, y_GRIEF_DE_NOSEL5))

y_GRIEF_DE_NOSEL =[x / 5 for x in grief_data]


#########################################

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp1/convergence.txt') as f:
    lines = f.readlines()

y_GRIEF_DE1 = [int(x) for x in lines]

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp2/convergence.txt') as f:
    lines = f.readlines()

y_GRIEF_DE2 = [int(x) for x in lines]

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp3/convergence.txt') as f:
    lines = f.readlines()

y_GRIEF_DE3 = [int(x) for x in lines]

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp4/convergence.txt') as f:
    lines = f.readlines()

y_GRIEF_DE4 = [int(x) for x in lines]

#with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp5/convergence.txt') as f:
#    lines = f.readlines()
#
#y_GRIEF_DE5 = [int(x) for x in lines]

grief_data = list( map(add, y_GRIEF_DE3, y_GRIEF_DE2))
grief_data = list( map(add, grief_data, y_GRIEF_DE3))
grief_data = list( map(add, grief_data, y_GRIEF_DE4))
#grief_data = list( map(add, grief_data, y_GRIEF_DE5))

y_GRIEF_DE =[x / 4 for x in grief_data]

##########################################




#with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp1/all_fits.txt') as f:
#    lines = f.readlines()
#
#count = 0
#data = []
#for line in lines:
#    result = re.sub(r'Gen [0-9]+', '', line)
#    result = re.sub(r'ind [0-9]+:', '', result)
#    if(result.replace("\n", "").split()):
#        data.append(result.replace("\n", "").split()[0])
#        data.append(result.replace("\n", "").split()[1])
#
#data= [int(x) for x in data]
#
#
#
#result = []
#for i in range(0,len(data),2):
#    result.append(max(data[:i+1]))
#
#
#with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp1/all_fits.txt') as f:
#    lines = f.readlines()
#
#count = 0
#data = []
#for line in lines:
#    r = re.sub(r'Gen [0-9]+', '', line)
#    r = re.sub(r'ind [0-9]+:', '', r)
#    if(r.replace("\n", "").split()):
#        data.append(r.replace("\n", "").split()[0])
#        data.append(r.replace("\n", "").split()[1])
#
#data= [int(x) for x in data]
#
#
#
#result = []
#for i in range(0,len(data),2):
#    result.append(max(data[:i+1]))

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp2/all_fits.txt') as f:
    lines = f.readlines()

count = 0
data = []
for line in lines:
    r = re.sub(r'Gen [0-9]+', '', line)
    r = re.sub(r'ind [0-9]+:', '', r)
    if(r.replace("\n", "").split()):
        data.append(r.replace("\n", "").split()[0])
        data.append(r.replace("\n", "").split()[1])

data= [int(x) for x in data]



result2 = []
for i in range(0,len(data),2):
    result2.append(max(data[:i+1]))

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp3/all_fits.txt') as f:
    lines = f.readlines()

count = 0
data = []
for line in lines:
    r = re.sub(r'Gen [0-9]+', '', line)
    r = re.sub(r'ind [0-9]+:', '', r)
    if(r.replace("\n", "").split()):
        data.append(r.replace("\n", "").split()[0])
        data.append(r.replace("\n", "").split()[1])

data= [int(x) for x in data]



result3 = []
for i in range(0,len(data),2):
    result3.append(max(data[:i+1]))

with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp4/all_fits.txt') as f:
    lines = f.readlines()

count = 0
data = []
for line in lines:
    r = re.sub(r'Gen [0-9]+', '', line)
    r = re.sub(r'ind [0-9]+:', '', r)
    if(r.replace("\n", "").split()):
        data.append(r.replace("\n", "").split()[0])
        data.append(r.replace("\n", "").split()[1])

data= [int(x) for x in data]



result4 = []
for i in range(0,len(data),2):
    result4.append(max(data[:i+1]))

#with open('/home/adriel/teste/GRIEF-DE (copy)/results/michigan/exp5/all_fits.txt') as f:
#    lines = f.readlines()
#
#count = 0
#data = []
#for line in lines:
#    r = re.sub(r'Gen [0-9]+', '', line)
#    r = re.sub(r'ind [0-9]+:', '', r)
#    if(r.replace("\n", "").split()):
#        data.append(r.replace("\n", "").split()[0])
#        data.append(r.replace("\n", "").split()[1])
#
#data= [int(x) for x in data]
#
#
#
#result5 = []
#for i in range(0,len(data),2):
#    result5.append(max(data[:i+1]))





#grief_data = list( map(add, result, result2))
#grief_data = list( map(add, grief_data, result3))
#grief_data = list( map(add, grief_data, result4))
#grief_data = list( map(add, grief_data, result5))

#y_GRIEF_DE_ =[x / 4 for x in grief_data]

#y_GRIEF_DE_ = y_GRIEF_DE_[:2000]

#plt.plot(y_GRIEF_DE)
#plt.show()
#
#print(len(result))
#exit()

#with open('./results/michigan/exp1/all_fits.txt') as f:
#    lines = f.readlines()

#count = 0
#data = []
#
#for line in lines:
#    r = re.sub(r'Gen [0-9]+', '', line)
#    r = re.sub(r'ind [0-9]+:', '', r)
#    if(r.replace("\n", "").split()):
#        data.append(r.replace("\n", "").split()[0])
#        data.append(r.replace("\n", "").split()[1])
#data= [int(x) for x in data]

#
#
#result = []
#for i in range(0,len(data),2):
#
#    result.append(max(data[:i+1]))
#with open('./results/michigan/exp2/all_fits.txt') as f:
#    lines = f.readlines()
#
#count = 0
#data = []
#
#for line in lines:
#    r = re.sub(r'Gen [0-9]+', '', line)
#    r = re.sub(r'ind [0-9]+:', '', r)
#    if(r.replace("\n", "").split()):
#        data.append(r.replace("\n", "").split()[0])
#        data.append(r.replace("\n", "").split()[1])
#data= [int(x) for x in data]

#
#
#result2 = []
#for i in range(0,len(data),2):
#
#    result2.append(max(data[:i+1]))
#with open('./results/michigan/exp3/all_fits.txt') as f:
#    lines = f.readlines()
#
#count = 0
#data = []
#
#for line in lines:
#    r = re.sub(r'Gen [0-9]+', '', line)
#    r = re.sub(r'ind [0-9]+:', '', r)
#    if(r.replace("\n", "").split()):
#        data.append(r.replace("\n", "").split()[0])
#        data.append(r.replace("\n", "").split()[1])
#data= [int(x) for x in data]

#
#
#result3 = []
#for i in range(0,len(data),2):
#
#    result3.append(max(data[:i+1]))
#with open('./results/michigan/exp4/all_fits.txt') as f:
#    lines = f.readlines()
#
#count = 0
#data = []
#
#for line in lines:
#    r = re.sub(r'Gen [0-9]+', '', line)
#    r = re.sub(r'ind [0-9]+:', '', r)
#    if(r.replace("\n", "").split()):
#        data.append(r.replace("\n", "").split()[0])
#        data.append(r.replace("\n", "").split()[1])
#data= [int(x) for x in data]

#
#
#result4 = []
#for i in range(0,len(data),2):
#
#    result4.append(max(data[:i+1]))
#
#grief_data = list( map(add, result, result2))
#
#grief_data = list( map(add, grief_data, result3))
#grief_data = list( map(add, grief_data, result4))
##grief_data = list( map(add, grief_data, result5))
#
#y_GRIEF_DE_30 =[x / 4 for x in grief_data]
#
#y_GRIEF_DE_30 = y_GRIEF_DE_30[:2000]



list_files = os.listdir("/home/adriel/repos/GRIEF/grief_history1")

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


list_files = os.listdir("/home/adriel/repos/GRIEF/grief_history2")

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

y_GRIEF2 = [x for _, x in sorted(zip(y_order, aux))]

list_files = os.listdir("/home/adriel/repos/GRIEF/grief_history3")

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

y_GRIEF3 = [x for _, x in sorted(zip(y_order, aux))]

list_files = os.listdir("/home/adriel/repos/GRIEF2/grief_history1")

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

y_GRIEF4 = [x for _, x in sorted(zip(y_order, aux))]

list_files = os.listdir("/home/adriel/repos/GRIEF3/grief_history1")

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

y_GRIEF5 = [x for _, x in sorted(zip(y_order, aux))]

grief_data = list( map(add, y_GRIEF, y_GRIEF2))
grief_data = list( map(add, grief_data, y_GRIEF3))
grief_data = list( map(add, grief_data, y_GRIEF4))
grief_data = list( map(add, grief_data, y_GRIEF5))

grief_data =[x / 5 for x in grief_data]

plt.title("Comparação do GRIEF e GRIEF-DE")
plt.plot(grief_data, label='GRIEF')
#plt.plot(y_GRIEF_DE_NOSEL, label = 'GRIEF-DE sem sel')
#plt.plot(y_GRIEF_DE_, label = 'GRIEF-DE sem sel (30 piores genes)')
plt.plot(y_GRIEF_DE_NOSEL, label = 'GRIEF-DE sem sel (10 piores genes)')
plt.legend()
plt.show()