#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np 
from math import sqrt
from scipy.integrate import solve_ivp

from os import listdir
from os.path import isfile, join


# In[2]:


def full_path(data_name, n):
    return '/home/dasha/picador/1.0rc1/bin/m=1/n5a25/BasicOutput/data/{:s}/{:06d}.bin'.format(data_name, n)


def load2DfromFile(dataname, f, ny):
    fd = open(full_path(dataname, f), 'rb')
    a = np.frombuffer(fd.read(), dtype=np.single)
    fd.close()
    return np.reshape(a, (ny, -1))

def files_in_dir(mypath):
    files = [file for file in listdir(mypath) if isfile(join(mypath,file))]
    return files

t = 20
m_e = 9.1/10**31
e = 1.6/10**19
c = 3*10**8 
m_i = m_e*1836

x = np.linspace(-2,10,600)
#                 ,data.shape[1]) #cm
p = np.linspace(-3.40581e-14,3.40581e-14,200) #г*см/с

data3 = load2DfromFile('electronphase',t,200)
data3 = load2DfromFile('ionphase',t,200)


x,y = np.meshgrid(x,p/m_i/c*1e-5)
ax = plt.contourf(x,y,data3,
                 cmap = 'magma',
                 levels = np.arange(0.005e16,2.5e16,0.1e16))
plt.colorbar(ax)
# plt.xlim(-1,10)
plt.ylim(-0.3,0.5)

plt.xlabel("$x$, $\mu$m")
plt.ylabel("$p_x/m_ic$")

# plt.savefig('fp_t20n5a25.pdf')


# In[9]:


def load2DfromFile(dataname, f, ny):
    fd = open(full_path(dataname, f), 'rb')
    a = np.frombuffer(fd.read(), dtype=np.single)
    fd.close()
    return(a)

t = 35
data = load2DfromFile('Ey',t,200)
x = np.linspace(-2,10,600)
p = np.linspace(-3.40581e-14,3.40581e-14,200) 
data = np.array(data)

plt.plot(x,data,linewidth=0.7)

data1 = load2DfromFile('Ex',t,200)
plt.plot(x,data1,linewidth=0.7)

n0 = 20
matrixsize_x = 600
x_max = 10/10**4
x_min = -2/10**4
Delta_x = (x_max-x_min)/matrixsize_x
n_cr = 1.1e21

ne = load2DfromFile('Electron1D',t,200)
n = ne/Delta_x/n_cr/n0
# plt.plot(x,n,linewidth=0.5,c='magenta')
# plt.ylim(0,55)
ni = load2DfromFile('Ion1D',t,200)
ni = ni/Delta_x/n_cr/n0
# plt.plot(x,ni,linewidth=0.5,c='b')
plt.xlim(-1,10)
plt.xlabel("$x$, $\mu$m")
plt.legend(("$E_y$","$E_x$"))
plt.savefig('pic_t35n5a25_.pdf')
# print(len(data))

