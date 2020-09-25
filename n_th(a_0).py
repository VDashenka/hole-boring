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


m_e = 9.1/10**31
e = 1.6/10**19
c = 3*10**8 
m_i = m_e*1836


# In[298]:


# распределение по импульсам
def full_path(data_name, n):
    return '/home/dasha/picador/1.0rc1/bin/n6a20m2.5/BasicOutput/data/{:s}/{:06d}.bin'.format(data_name, n)


def load2DfromFile(dataname, f, ny):
    fd = open(full_path(dataname, f), 'rb')
    a = np.frombuffer(fd.read(), dtype=np.single)
    fd.close()
    return np.reshape(a, (ny, -1))

def files_in_dir(mypath):
    files = [file for file in listdir(mypath) if isfile(join(mypath,file))]
    return files

data = load2DfromFile('ionphase',100,200)
data = np.array(data)
print(data.shape)

x = np.linspace(-0.0002,0.001,data.shape[1]) #cм
p = np.linspace(-1.48021e-14,1.48021e-14,data.shape[0]) #г*см/с

f = np.trapz(data,x,axis=1)
print(f.shape)
print(len(x),len(p))

plt.plot(p,f)
p_i = p[np.argmax(f)]/m_i/c*1e-5
print(p_i)


# In[4]:


# энергия теория
a_list = [10,15,20,25]
n_list = [4,6,8,10]

P_i = np.zeros((len(n_list),len(a_list)))
W_i = np.zeros((len(n_list),len(a_list)))
V_i = np.zeros((len(n_list),len(a_list)))


for i in range(len(a_list)):
    n_0 = n_list[i] 
    for j in range(len(n_list)):
        
        m = 1.0/(1836) #\mu
        a_pad = a_list[j]

        print(n_0,'n_0 задано')
        print(a_pad,'a_pad')
        
        b = a_pad/(sqrt(n_0/m)+a_pad)
        ge = 1/sqrt(1-b**2)
        gi = ge
        print(ge,'gamma')

        p = n_0*b

        #начальное значение
        a_0=0.001

        #граничные точки
        z_0=0
        z_f=2.5

        l=p/sqrt(ge**2-1)

        phi_0=0
        a1_0=sqrt(l-1)*a_0
        phi1_0=2*sqrt(l-1)*phi_0


        yy = np.ndarray(shape=(4,1))

        def f(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        p*((ge+y1)/np.sqrt(((y1+ge)**2-(y3)**2)-1)-(gi-y1*m)/np.sqrt((gi-y1*m)**2-(m*y3)**2-1)),
                        y4,
                        -y3+p*y3/np.sqrt(((y1+ge)**2-(y3)**2)-1) ]

        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        k= 0
        try:
            k=solve_ivp(f,(z_0,z_f),y_0)#, max_step=1e-4)
        except:
            print("error")

        phi_p = k.y[0]
        E_p=k.y[1]
        a_p = k.y[2]
        z_p=k.t



        #ионная область
        z_0i=k.t[-1]
        z_fi=z_f+0.5
#         print(z_0i,'граница ионы-плазма')

        phi_0=k.y[0][-1]
        phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
        a_0=k.y[2][-1]
        a1_0=sqrt(4*a_pad**2-a_0**2)
        a1_0_1=k.y[3][-1]
        

        def f_i(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                    y4,
                    -y3 ]


        y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
        except:
            print("error")


        phi_i = ki.y[0]
        E_i=ki.y[1]
        # print(E_i)
        a_i = ki.y[2]
        z_i=ki.t


        #вакуумная область
        z_0v=ki.t[-1]
        z_fv=z_fi+1.5

        phi_0=ki.y[0][-1]
        phi1_0=ki.y[1][-1]
        a_0=ki.y[2][-1]
#         print(a_0,'а_0 на гр вак-ион  1')
        a1_0=ki.y[3][-1]

        def f_v(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    0,
                    y4,
                    -y3 ]


        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
        except:
            print("error")

        phi_v = kv.y[0]
        E_v=kv.y[1]
        a_v = kv.y[2]
        z_v=kv.t
        
        a=np.hstack((a_p,a_i,a_v))
        z=np.hstack((z_p,z_i,z_v))
        #plt.plot(z,abs(a))


        E=np.hstack((E_p,E_i,E_v))
        phi=np.hstack((phi_p,phi_i,phi_v))
        #plt.plot(z,E)
        #plt.plot(z,phi)
        
        
        if E_i[-1]<=0:
            print('ddd')
            a=0.1 #точность
            for i in range(0,len(E_i)):
                if E_i[i]>=0:
                    d=i
                # print(d,z[d],E[d])
                    if E_i[i]<E_i[d]:
                        d=i
       
        
        #ионная область
            z_0i=k.t[-1]
            z_fi=z_i[d]
        # print(z_0i,'граница ионы-плазма')

            phi_0=k.y[0][-1]
            phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
            a_0=k.y[2][-1]
            a1_0=sqrt(4*a_pad**2-a_0**2)
            a1_0_1=k.y[3][-1]


            def f_i(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                        y4,
                        -y3 ]


            y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
            except:
                print("error")


            phi_i = ki.y[0]
            E_i=ki.y[1]
            a_i = ki.y[2]
            z_i=ki.t


        #вакуумная область
            z_0v=ki.t[-1]
            z_fv=z_fi+1.5

            phi_0=ki.y[0][-1]
            phi1_0=ki.y[1][-1]
            a_0=ki.y[2][-1]
      
            a1_0=ki.y[3][-1]

            def f_v(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        0,
                        y4,
                        -y3 ]


            y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
            except:
                print("error")

            phi_v = kv.y[0]
            E_v=kv.y[1]
            a_v = kv.y[2]
            z_v=kv.t
        
            a=np.hstack((a_p,a_i,a_v))
            z=np.hstack((z_p,z_i,z_v))
            #plt.plot(z,abs(a))

            E=np.hstack((E_p,E_i,E_v))
            phi=np.hstack((phi_p,phi_i,phi_v))
        
        m_i = m_e/m
        w_i = 2*m_i*c**2*ge**2*b**2
        w_i = w_i/e/10**6 #MeV
        print(w_i,'энергия, МэВ')
        v_i = 2*b*c/(1+b**2)
        print(v_i,'скорость ионов')
        
        W_i[i][j] = w_i
        V_i[i][j] = v_i


# In[278]:


print(W_i)


# In[281]:


x,y = np.meshgrid(a_list,n_list)

ax1 = plt.contourf(x,y,W_i,cmap='inferno')
#                 levels = np.arange(0,180,2))
plt.colorbar(ax1)
plt.xlabel('$a_0$')
plt.ylabel('$n_0$')


# In[5]:


# n4a10 (-2.34042e-14,2.34042e14)
# n4a15 (-3.51063e-14,3.51063e-14)
# n4a20 (-4.68084e-14,4.68084e-14)
# n4a25 (-5.85106e-14,5.85106e-14)

# n5a10 -2.09334e-14
# n5a15 -3.14001e-14
# n5a20 -4.18667e-14
# n5a25 -5.23334e-14
#        [-2.09334e-14,-3.14001e-14,-4.18667e-14,-5.23334e-14],


# n6a10 (-1.91095e-14,1.91095e-14)
# n6a15 (-2.86642e-14,2.86642e-14)
# n6a20 (-3.82189e-14,3.82189e-14)
# n6a25 (-4.77737e-14,4.77737e-14)

# n8a10 (-1.65493e-14,1.65493e-14)
# n8a15 (-2.48239e-14,2.48239e-14)
# n8a20 (-3.30986e-14,3.30986e-14)
# n8a25 (-4.13732e-14,4.13732e-14)

# n10a10 (-1.48021e-14,1.48021e-14)
# n10a15 (-2.22032e-14,2.22032e-14)
# n10a20 (-2.96043e-14,2.96043e-14)
# n10a25 (-3.70053e-14,3.70053e-14)

dp = [[-2.34042e-14,-1.91095e-14,-1.65493e-14,-1.48021e-14],
    [-3.51063e-14,-2.86642e-14,-2.48239e-14,-2.22032e-14],
    [-4.68084e-14,-3.82189e-14,-3.30986e-14,-2.96043e-14],
    [-5.85106e-14,-4.77737e-14,-4.13732e-14,-3.70053e-14]]

dp1 = [[-2.34042e-14,-3.51063e-14,-4.68084e-14,-5.85106e-14],
      [-1.91095e-14,-2.86642e-14,-3.82189e-14,-4.77737e-14],
      [-1.65493e-14,-2.48239e-14,-3.30986e-14,-4.13732e-14],
      [-1.48021e-14,-2.22032e-14,-2.96043e-14,-3.70053e-14]]


# In[6]:


n_0 = n_list
a_0 = a_list

x = np.linspace(-0.0002,0.001,600)
# p = np.linspace(-3.40581e-14,3.40581e-14,200)

P = np.zeros((len(n_0),len(a_0)))
Pa = np.zeros((len(n_0),len(a_0)))

V_pic = np.zeros((len(n_0),len(a_0)))

E = np.zeros((len(n_0),len(a_0)))
delta = np.zeros((len(n_0),len(a_0)))
delta1 = np.zeros((len(n_0),len(a_0)))


df = ['n4','n6','n8','n10']

for i in range(len(n_0)):
    fail = files_in_dir('/home/dasha/picador/1.0rc1/bin/t_100/'+df[i])
    fail.sort()
    print(fail,df[i])
    for j in range(len(a_0)):
        fd = open('/home/dasha/picador/1.0rc1/bin/t_100/'+df[i]+'/'+fail[j], 'rb')
        print(fail[j])
        
        a = np.frombuffer(fd.read(), dtype=np.single)
        fd.close()
        data = np.reshape(a, (200, -1))
        
        p = np.linspace(dp1[i][j],-dp1[i][j],200)
        print(dp1[i][j])

        data = np.array(data)
        f = np.trapz(data,x,axis=1)
        f.shape
        plt.plot(p,f)
        plt.show()
        
        p_i = p[np.argmax(f)]/m_i/c*1e-5
        print(p_i,'импульс')

        gamma = sqrt(p_i**2+1)
        v_i_pic = c*sqrt(1-1/gamma/gamma)
        V_pic[i][j] = v_i_pic
        print(v_i_pic/c,'скорость/c')
        
        E_i = m_i*c**2*(gamma-1)/e/10**6 #MeV
        print(E_i,'энергия')
        E[i][j] = E_i
        delta[i][j] = E[i][j]/W_i[i][j] #pic/th

        P[i][j] = p_i        

print(E)


# In[284]:


print(E)
print(W_i)
print(delta)


# In[271]:


print(E[0][1])
print(W_i[1][0])


# In[7]:


a_list = [10,15,20,25]
n_list = [4,6,8,10]
n_0 = n_list
a_0 = a_list

print(len(n_0),len(n_0),delta.shape)
x,y = np.meshgrid(a_0,n_0)

ax4 = plt.contourf(x,y,delta,levels = 100,
#                    cmap='spring',
                   cmap='autumn')
#                    cmap='RdYlBu',#pic
#                    np.arange(0,1.45,0.05))
plt.colorbar(ax4)
plt.xlabel('$a_0$')
plt.ylabel('$n_0$')

plt.contour(x, y, delta,
                levels = [0.998],
                colors = 'navy')
# plt.title('$W_{pic}/W_{th}$')
# plt.savefig('n_th5.pdf')


# In[64]:


print(a_0,n_0)
# a_0 = [10,15,20,25]
# n_0 = [4,6,8,10]
a_list = [10,15,20,25]
n_list = [4,6,8,10]
n_0 = n_list
a_0 = a_list

x,y = np.meshgrid(n_0,a_0)

ax4 = plt.contourf(x,y,delta,levels = 100)
#                    cmap='spring',
#                    cmap='autumn',
#                    cmap='RdYlBu',#pic
#                    np.arange(0,1.45,0.05))
plt.colorbar(ax4)
plt.xlabel('$n_0$')
plt.ylabel('$a_0$')

plt.contour(x, y, delta,
                levels = [0.998],
                colors = 'navy')
plt.title('$W_{pic}/W_{th}$')

plt.savefig('WW.pdf')
# n = [10,9,8,]
# a = [20,15,11]
# plt.plot(n,a)


# In[151]:


# энергия теория
a_list = [10,15,20,25]
n_list = [4,6,8,10]

P_i5 = np.zeros((len(n_list),len(a_list)))
W_i5 = np.zeros((len(n_list),len(a_list)))
V_i5 = np.zeros((len(n_list),len(a_list)))


for ii in range(len(n_list)):
    n_0 = n_list[ii] 
    for j in range(len(a_list)):
        
        m = 1.0/(1836)/5 #\mu
        a_pad = a_list[j]

        print(n_0,'n_0 задано')
        print(a_pad,'a_pad')
        
        b = a_pad/(sqrt(n_0/m)+a_pad)
        ge = 1/sqrt(1-b**2)
        gi = ge
        print(ge,'gamma')

        p = n_0*b

        #начальное значение
        a_0=0.001

        #граничные точки
        z_0=0
        z_f=2.5

        l=p/sqrt(ge**2-1)

        phi_0=0
        a1_0=sqrt(l-1)*a_0
        phi1_0=2*sqrt(l-1)*phi_0


        yy = np.ndarray(shape=(4,1))

        def f(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        p*((ge+y1)/np.sqrt(((y1+ge)**2-(y3)**2)-1)-(gi-y1*m)/np.sqrt((gi-y1*m)**2-(m*y3)**2-1)),
                        y4,
                        -y3+p*y3/np.sqrt(((y1+ge)**2-(y3)**2)-1) ]

        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        k= 0
        try:
            k=solve_ivp(f,(z_0,z_f),y_0)#, max_step=1e-4)
        except:
            print("error")

        phi_p = k.y[0]
        E_p=k.y[1]
        a_p = k.y[2]
        z_p=k.t



        #ионная область
        z_0i=k.t[-1]
        z_fi=z_f+0.5
#         print(z_0i,'граница ионы-плазма')

        phi_0=k.y[0][-1]
        phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
        a_0=k.y[2][-1]
        a1_0=sqrt(4*a_pad**2-a_0**2)
        a1_0_1=k.y[3][-1]
        

        def f_i(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                    y4,
                    -y3 ]


        y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
        except:
            print("error")


        phi_i = ki.y[0]
        E_i=ki.y[1]
        # print(E_i)
        a_i = ki.y[2]
        z_i=ki.t


        #вакуумная область
        z_0v=ki.t[-1]
        z_fv=z_fi+1.5

        phi_0=ki.y[0][-1]
        phi1_0=ki.y[1][-1]
        a_0=ki.y[2][-1]
#         print(a_0,'а_0 на гр вак-ион  1')
        a1_0=ki.y[3][-1]

        def f_v(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    0,
                    y4,
                    -y3 ]


        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
        except:
            print("error")

        phi_v = kv.y[0]
        E_v=kv.y[1]
        a_v = kv.y[2]
        z_v=kv.t
        
        a=np.hstack((a_p,a_i,a_v))
        z=np.hstack((z_p,z_i,z_v))
        #plt.plot(z,abs(a))


        E=np.hstack((E_p,E_i,E_v))
        phi=np.hstack((phi_p,phi_i,phi_v))
        #plt.plot(z,E)
        #plt.plot(z,phi)
        
        
        if E_i[-1]<=0:
            print('ddd')
            a=0.1 #точность
            for i in range(0,len(E_i)):
                if E_i[i]>=0:
                    d=i
                # print(d,z[d],E[d])
                    if E_i[i]<E_i[d]:
                        d=i
       
        
        #ионная область
            z_0i=k.t[-1]
            z_fi=z_i[d]
        # print(z_0i,'граница ионы-плазма')

            phi_0=k.y[0][-1]
            phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
            a_0=k.y[2][-1]
            a1_0=sqrt(4*a_pad**2-a_0**2)
            a1_0_1=k.y[3][-1]


            def f_i(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                        y4,
                        -y3 ]


            y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
            except:
                print("error")


            phi_i = ki.y[0]
            E_i=ki.y[1]
            a_i = ki.y[2]
            z_i=ki.t


        #вакуумная область
            z_0v=ki.t[-1]
            z_fv=z_fi+1.5

            phi_0=ki.y[0][-1]
            phi1_0=ki.y[1][-1]
            a_0=ki.y[2][-1]
      
            a1_0=ki.y[3][-1]

            def f_v(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        0,
                        y4,
                        -y3 ]


            y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
            except:
                print("error")

            phi_v = kv.y[0]
            E_v=kv.y[1]
            a_v = kv.y[2]
            z_v=kv.t
        
            a=np.hstack((a_p,a_i,a_v))
            z=np.hstack((z_p,z_i,z_v))
            #plt.plot(z,abs(a))

            E=np.hstack((E_p,E_i,E_v))
            phi=np.hstack((phi_p,phi_i,phi_v))
        
        m_i = m_e/m
        w_i = 2*m_i*c**2*ge**2*b**2
        w_i = w_i/e/10**6 #MeV
        print(w_i,'энергия, МэВ')
        v_i = 2*b*c/(1+b**2)
        print(v_i,'скорость ионов')
        
        W_i5[ii][j] = w_i
        V_i5[ii][j] = v_i


# In[152]:


x,y = np.meshgrid(a_list,n_list)

ax1 = plt.contourf(x,y,W_i5,cmap='inferno')
#                 levels = np.arange(0,180,2))
plt.colorbar(ax1)
plt.xlabel('$a_0$')
plt.ylabel('$n_0$')


# In[119]:


# n4a10 -3.51063e-14

dp5 = [[-3.51063e-14,-3.51063e-14,-4.68084e-14,-5.85106e-14],
      [-2.86642e-14,-4.29963e-14,-5.73284e-14,-7.16605e-14],
      [-2.48239e-14,-3.72359e-14,-4.96478e-14,-6.20598e-14],
      [-2.22032e-14,-3.33048e-14,-4.44064e-14,-5.5508e-14]]


# In[159]:


n_0 = n_list
a_0 = a_list

x = np.linspace(-0.0002,0.001,600)

P5 = np.zeros((len(n_0),len(a_0)))
Pa5 = np.zeros((len(n_0),len(a_0)))

V_pic5 = np.zeros((len(n_0),len(a_0)))

E5 = np.zeros((len(n_0),len(a_0)))
delta5 = np.zeros((len(n_0),len(a_0)))

df = ['n4','n6','n8','n10']

for i in range(len(n_0)):
    fail = files_in_dir('/home/dasha/picador/1.0rc1/bin/m5/'+df[i])
    fail.sort()
    print(fail,df[i])
    for j in range(len(a_0)):
        fd = open('/home/dasha/picador/1.0rc1/bin/m5/'+df[i]+'/'+fail[j], 'rb')
        print(fail[j])
        
        a = np.frombuffer(fd.read(), dtype=np.single)
        fd.close()
        data = np.reshape(a, (200, -1))
        
        p = np.linspace(dp5[i][j],-dp5[i][j],200)
        print(dp5[i][j])

        data = np.array(data)
        f = np.trapz(data,x,axis=1)
        f.shape
        p = p[120:len(p)]
        f = f[120:len(f)]
#         plt.plot(p,f)
#         plt.show()
#         print(p)
#         print(p)
#         if df[i] == 'n4' and fail[j] == 'a10': 
# #         or (df[i] == 'n6' and fail[j] == 'a10'):
#             p = p[120:len(p)]
#             f = f[120:len(f)]
#             plt.plot(p,f)
#             plt.show()

#             p_i = p[np.argmax(f)]/m_i/c*1e-5
#             print(p_i,'импульс')

#             gamma = sqrt(p_i**2+1)
#             v_i_pic = c*sqrt(1-1/gamma/gamma)
#             V_pic[i][j] = v_i_pic
#             print(v_i_pic/c,'скорость/c')
        
#             E_i = m_i*c**2*(gamma-1)/e/10**6 #MeV
#             print(E_i,'энергия')
#             E5[i][j] = E_i
#             delta5[i][j] = E5[i][j]/W_i5[i][j] #pic/th

#             P5[i][j] = p_i   
        
#             p_i = p[np.argmax(f)]/m_i/c*1e-5
#             print(p_i,'импульс')

#             gamma = sqrt(p_i**2+1)
#             v_i_pic = c*sqrt(1-1/gamma/gamma)
#             V_pic[i][j] = v_i_pic
#             print(v_i_pic/c,'скорость/c')
        
#             E_i = m_i*c**2*(gamma-1)/e/10**6 #MeV
#             print(E_i,'энергия')
#             E5[i][j] = E_i
#             delta5[i][j] = E5[i][j]/W_i5[i][j] #pic/th

#             P5[i][j] = p_i    
            
#         else:
        p_i = p[np.argmax(f)]/m_i/c*1e-5
        print(p_i,'импульс')
        plt.plot(p,f)


        gamma = sqrt(p_i**2+1)
        v_i_pic = c*sqrt(1-1/gamma/gamma)
        V_pic[i][j] = v_i_pic
        print(v_i_pic/c,'скорость/c')
        
        E_i = m_i*c**2*(gamma-1)/e/10**6 #MeV
        print(E_i,'энергия')
        E5[i][j] = E_i
        delta5[i][j] = E5[i][j]/W_i5[i][j] #pic/th

        P5[i][j] = p_i   

print(E5)


# In[153]:


print(E5)
print(W_i5)


# In[261]:


# энергия теория
a_list = [10,15,20,25]
n_list = [4,6,8,10]

P_i25 = np.zeros((len(n_list),len(a_list)))
W_i25 = np.zeros((len(n_list),len(a_list)))
V_i25 = np.zeros((len(n_list),len(a_list)))


for ii in range(len(n_list)):
    n_0 = n_list[ii] 
    for j in range(len(a_list)):
        
        m = 1.0/(1836)/2.5 #\mu
        a_pad = a_list[j]

        print(n_0,'n_0 задано')
        print(a_pad,'a_pad')
        
        b = a_pad/(sqrt(n_0/m)+a_pad)
        ge = 1/sqrt(1-b**2)
        gi = ge
        print(ge,'gamma')

        p = n_0*b

        #начальное значение
        a_0=0.001

        #граничные точки
        z_0=0
        z_f=2.5

        l=p/sqrt(ge**2-1)

        phi_0=0
        a1_0=sqrt(l-1)*a_0
        phi1_0=2*sqrt(l-1)*phi_0


        yy = np.ndarray(shape=(4,1))

        def f(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        p*((ge+y1)/np.sqrt(((y1+ge)**2-(y3)**2)-1)-(gi-y1*m)/np.sqrt((gi-y1*m)**2-(m*y3)**2-1)),
                        y4,
                        -y3+p*y3/np.sqrt(((y1+ge)**2-(y3)**2)-1) ]

        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        k= 0
        try:
            k=solve_ivp(f,(z_0,z_f),y_0)#, max_step=1e-4)
        except:
            print("error")

        phi_p = k.y[0]
        E_p=k.y[1]
        a_p = k.y[2]
        z_p=k.t



        #ионная область
        z_0i=k.t[-1]
        z_fi=z_f+0.5
#         print(z_0i,'граница ионы-плазма')

        phi_0=k.y[0][-1]
        phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
        a_0=k.y[2][-1]
        a1_0=sqrt(4*a_pad**2-a_0**2)
        a1_0_1=k.y[3][-1]
        

        def f_i(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                    y4,
                    -y3 ]


        y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
        except:
            print("error")


        phi_i = ki.y[0]
        E_i=ki.y[1]
        # print(E_i)
        a_i = ki.y[2]
        z_i=ki.t


        #вакуумная область
        z_0v=ki.t[-1]
        z_fv=z_fi+1.5

        phi_0=ki.y[0][-1]
        phi1_0=ki.y[1][-1]
        a_0=ki.y[2][-1]
#         print(a_0,'а_0 на гр вак-ион  1')
        a1_0=ki.y[3][-1]

        def f_v(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    0,
                    y4,
                    -y3 ]


        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
        except:
            print("error")

        phi_v = kv.y[0]
        E_v=kv.y[1]
        a_v = kv.y[2]
        z_v=kv.t
        
        a=np.hstack((a_p,a_i,a_v))
        z=np.hstack((z_p,z_i,z_v))
        #plt.plot(z,abs(a))


        E=np.hstack((E_p,E_i,E_v))
        phi=np.hstack((phi_p,phi_i,phi_v))
        #plt.plot(z,E)
        #plt.plot(z,phi)
        
        
        if E_i[-1]<=0:
            print('ddd')
            a=0.1 #точность
            for i in range(0,len(E_i)):
                if E_i[i]>=0:
                    d=i
                # print(d,z[d],E[d])
                    if E_i[i]<E_i[d]:
                        d=i
       
        
        #ионная область
            z_0i=k.t[-1]
            z_fi=z_i[d]
        # print(z_0i,'граница ионы-плазма')

            phi_0=k.y[0][-1]
            phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
            a_0=k.y[2][-1]
            a1_0=sqrt(4*a_pad**2-a_0**2)
            a1_0_1=k.y[3][-1]


            def f_i(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                        y4,
                        -y3 ]


            y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
            except:
                print("error")


            phi_i = ki.y[0]
            E_i=ki.y[1]
            a_i = ki.y[2]
            z_i=ki.t


        #вакуумная область
            z_0v=ki.t[-1]
            z_fv=z_fi+1.5

            phi_0=ki.y[0][-1]
            phi1_0=ki.y[1][-1]
            a_0=ki.y[2][-1]
      
            a1_0=ki.y[3][-1]

            def f_v(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        0,
                        y4,
                        -y3 ]


            y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
            except:
                print("error")

            phi_v = kv.y[0]
            E_v=kv.y[1]
            a_v = kv.y[2]
            z_v=kv.t
        
            a=np.hstack((a_p,a_i,a_v))
            z=np.hstack((z_p,z_i,z_v))
            #plt.plot(z,abs(a))

            E=np.hstack((E_p,E_i,E_v))
            phi=np.hstack((phi_p,phi_i,phi_v))
        
        m_i = m_e/m
        w_i = 2*m_i*c**2*ge**2*b**2
        w_i = w_i/e/10**6 #MeV
        print(w_i,'энергия, МэВ')
        v_i = 2*b*c/(1+b**2)
        print(v_i,'скорость ионов')
        
        W_i25[ii][j] = w_i
        V_i25[ii][j] = v_i


# In[262]:


x,y = np.meshgrid(a_list,n_list)

ax1 = plt.contourf(x,y,W_i25,cmap='inferno')
#                 levels = np.arange(0,180,2))
plt.colorbar(ax1)
plt.xlabel('$a_0$')
plt.ylabel('$n_0$')


# In[260]:


dp25 = [[-2.92553e-14,-4.38829e-14,-5.85106e-14,-7.31382e-14],
       [-2.38868e-14,-3.58302e-14,-4.77737e-14,-5.97171e-14],
       [-2.06866e-14,-3.10299e-14,-4.13732e-14,-5.17165e-14],
       [-1.85027e-14,-2.7754e-14,-3.70053e-14,-4.62567e-14]]


# In[337]:


n_0 = n_list
a_0 = a_list

x = np.linspace(-0.0002,0.001,600)

P25 = np.zeros((len(n_0),len(a_0)))
Pa25 = np.zeros((len(n_0),len(a_0)))

V_pic25 = np.zeros((len(n_0),len(a_0)))

E25 = np.zeros((len(n_0),len(a_0)))
delta25 = np.zeros((len(n_0),len(a_0)))

df = ['n4','n6','n8','n10']

for i in range(len(n_0)):
    fail = files_in_dir('/home/dasha/picador/1.0rc1/bin/m2.5/'+df[i])
    fail.sort()
    print(fail,df[i])
    for j in range(len(a_0)):
        fd = open('/home/dasha/picador/1.0rc1/bin/m2.5/'+df[i]+'/'+fail[j], 'rb')
        print(fail[j])
        
        a = np.frombuffer(fd.read(), dtype=np.single)
        fd.close()
        data = np.reshape(a, (200, -1))
        
        p = np.linspace(dp5[i][j],-dp5[i][j],200)
        print(dp5[i][j])

        data = np.array(data)
        f = np.trapz(data,x,axis=1)
        f.shape
        p = p[120:len(p)]
        f = f[120:len(f)]

        p_i = p[np.argmax(f)]/m_i/c*1e-5
        print(p_i,'импульс')
        plt.plot(p,f)


        gamma = sqrt(p_i**2+1)
        v_i_pic = c*sqrt(1-1/gamma/gamma)
        V_pic25[i][j] = v_i_pic
        print(v_i_pic/c,'скорость/c')
        
        E_i = m_i*c**2*(gamma-1)/e/10**6 #MeV
        print(E_i,'энергия')
        E25[i][j] = E_i
        delta25[i][j] = E25[i][j]/W_i25[i][j] #pic/th

        P25[i][j] = p_i   

print(E25)


# In[300]:


x,y = np.meshgrid(a_0,n_0)

ax4 = plt.contourf(x,y,delta25,levels = 100,
#                    cmap='spring',
                   cmap='autumn')
#                    cmap='RdYlBu',#pic
#                    np.arange(0,1.45,0.05))
plt.colorbar(ax4)

plt.contour(x, y, delta25,
                levels = [1],
                colors = 'navy')


# In[354]:


a_list = [10,15,20,25]
n_list = [4,6,8,10]
n_0 = n_list
a_0 = a_list

print(len(n_0),len(n_0),delta.shape)
x,y = np.meshgrid(a_0,n_0)

# ax4 = plt.contourf(x,y,delta5,levels = 100,
#                    cmap='spring',
#                    cmap='autumn')
#                    cmap='RdYlBu',#pic
#                    np.arange(0,1.45,0.05))
# plt.colorbar(ax4)
plt.xlabel('$a_0$')
plt.ylabel('$n_0$')

plt.contour(x, y, delta5,
                levels = [0.9],
                colors = 'green')

plt.contour(x, y, delta25,
                levels = [0.99],
                colors = 'navy')

plt.contour(x, y, delta,
                levels = [0.9],
                colors = 'magenta')

plt.text(16, 7.1, '$\mu=5$', withdash=False,fontdict={'family':'verdana'},size="x-large")
plt.text(19.5, 6.2, '$\mu=2.5$', withdash=False,fontdict={'family':'verdana'},size="x-large")
plt.text(22, 5, '$\mu=1$', withdash=False,fontdict={'family':'verdana'},size="x-large")


# plt.title('$W_{pic}/W_{th}$')
plt.savefig('dif_mu.pdf')


# In[8]:


# P_cr/a_0

A_pad_1 = [2,6,10,14,18,22,25]
P_cr = [] 
n_list_1 = [2,5,8,11,14,17]

P = np.zeros((len(n_list_1),len(A_pad_1)))
P_A = np.zeros((len(n_list_1),len(A_pad_1)))


for iii in range(len(n_list_1)):
    n_0 = n_list_1[iii] 

    print(a_pad,'a_pad задано')
    for j in range(len(A_pad_1)): 
        
        m = 1.0/(1836) #\mu
        a_pad = A_pad_1[j]

        print(n_0,'n_0 задано')
       
        b = a_pad/(sqrt(n_0/m)+a_pad)
        ge = 1/sqrt(1-b**2)
        gi = ge
        print(ge,'gamma')

        p = n_0*b

        #начальное значение
        a_0=0.001

        #граничные точки
        z_0=0
        z_f=2.5

        l=p/sqrt(ge**2-1)

        phi_0=0
        a1_0=sqrt(l-1)*a_0
        phi1_0=2*sqrt(l-1)*phi_0


        yy = np.ndarray(shape=(4,1))

        def f(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        p*((ge+y1)/np.sqrt(((y1+ge)**2-(y3)**2)-1)-(gi-y1*m)/np.sqrt((gi-y1*m)**2-(m*y3)**2-1)),
                        y4,
                        -y3+p*y3/np.sqrt(((y1+ge)**2-(y3)**2)-1) ]

        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        k= 0
        try:
            k=solve_ivp(f,(z_0,z_f),y_0)#, max_step=1e-4)
        except:
            print("error")

        phi_p = k.y[0]
        E_p=k.y[1]
        a_p = k.y[2]
        z_p=k.t



        #ионная область
        z_0i=k.t[-1]
        z_fi=z_f+0.5
        print(z_0i,'граница ионы-плазма')

        phi_0=k.y[0][-1]
        phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
        a_0=k.y[2][-1]
        a1_0=sqrt(4*a_pad**2-a_0**2)
        a1_0_1=k.y[3][-1]
        
        print(a_0,'а на гр ионы-пл')
        print(a_pad,'задан a_pad')

        def f_i(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                    y4,
                    -y3 ]


        y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
        except:
            print("error")


        phi_i = ki.y[0]
        E_i=ki.y[1]
        # print(E_i)
        a_i = ki.y[2]
        z_i=ki.t


        #вакуумная область
        z_0v=ki.t[-1]
        z_fv=z_fi+1.5

        phi_0=ki.y[0][-1]
        phi1_0=ki.y[1][-1]
        a_0=ki.y[2][-1]
        print(a_0,'а_0 на гр вак-ион  1')
        a1_0=ki.y[3][-1]

        def f_v(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    0,
                    y4,
                    -y3 ]


        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
        except:
            print("error")

        phi_v = kv.y[0]
        E_v=kv.y[1]
        a_v = kv.y[2]
        z_v=kv.t
        
        a=np.hstack((a_p,a_i,a_v))
        z=np.hstack((z_p,z_i,z_v))
        #plt.plot(z,abs(a))


        E=np.hstack((E_p,E_i,E_v))
        phi=np.hstack((phi_p,phi_i,phi_v))
        #plt.plot(z,E)
        #plt.plot(z,phi)
        
        
        if E_i[-1]<=0:
            print('ddd')
            a=0.1 #точность
            for i in range(0,len(E_i)):
                if E_i[i]>=0:
                    d=i
                # print(d,z[d],E[d])
                    if E_i[i]<E_i[d]:
                        d=i
       

        #print(d,z[d],E[d])
        
        #ионная область
            z_0i=k.t[-1]
            z_fi=z_i[d]
        # print(z_0i,'граница ионы-плазма')

            phi_0=k.y[0][-1]
            phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
            a_0=k.y[2][-1]
            a1_0=sqrt(4*a_pad**2-a_0**2)
            a1_0_1=k.y[3][-1]

        #print(a1_0,"а' из связи а и а'")
        #print(a1_0_1, "a' из построенного решения")


            def f_i(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                        y4,
                        -y3 ]


            y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
            except:
                print("error")


            phi_i = ki.y[0]
            E_i=ki.y[1]
            a_i = ki.y[2]
            z_i=ki.t


        #вакуумная область
            z_0v=ki.t[-1]
            z_fv=z_fi+1.5

            phi_0=ki.y[0][-1]
            phi1_0=ki.y[1][-1]
            a_0=ki.y[2][-1]
        #print(a_0)
            print(a_0,'а_0 на гр вак-ион  2')
        #a_pad=a_0
            a1_0=ki.y[3][-1]

            def f_v(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        0,
                        y4,
                        -y3 ]


            y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
            except:
                print("error")

            phi_v = kv.y[0]
            E_v=kv.y[1]
            a_v = kv.y[2]
            z_v=kv.t
        
            a=np.hstack((a_p,a_i,a_v))
            z=np.hstack((z_p,z_i,z_v))
            #plt.plot(z,abs(a))


            E=np.hstack((E_p,E_i,E_v))
            phi=np.hstack((phi_p,phi_i,phi_v))

        
        # ищем седло
        a=0.005 #точность
        z=np.hstack((z_i,z_v))
        for i in range(0,len(z_i)):
            if ((a_i[i]*ki.y[3][i])/np.sqrt(1+(a_i[i])**2))-E_i[i]<a and ((a_i[i]*ki.y[3][i])/np.sqrt(1+(a_i[i])**2))-E_i[i]>-a:
                print(i)
            
        for i in range(0,len(z_v)):
            if ((a_v[i]*kv.y[3][i])/np.sqrt(1+(a_v[i])**2))-E_v[i]<a and ((a_v[i]*kv.y[3][i])/np.sqrt(1+(a_v[i])**2))-E_v[i]>-a:
        #         print(i,"в вак")
                g=i
        #print(g)  
        
        s=g #индекс, соответствующий седловой точке
        h_cr=np.sqrt(1+(a_v[s])**2)-phi_v[s]
        p_cr=np.sqrt((h_cr+phi_i[0])**2-(a_i[0])**2-1)
        #print(p_cr,'критический импульс электрона на границе с ионами')
        
        P_cr.append(p_cr)
        P[iii][j] = p_cr
        P_A[iii][j] = p_cr/a_pad


# In[174]:


print(n_list_1)
# X,Y = np.meshgrid(A_pad,n_list)
# ax = plt.contourf(X,Y,P_A,levels = 100)
#                  cmap='autumn')
#                   np.arange(0,3,0.1))

# plt.colorbar(ax)
plt.contour(x, y, P_A_1,
                levels = [1],
                colors = 'cyan')

a_0 = [10,15,20,25]
n_0 = [4,6,8,10]

x,y = np.meshgrid(a_0,n_0)
# ax4 = plt.contourf(x,y,delta,levels = 100)
#                    cmap='spring',
#                    cmap='autumn',
#                    cmap='RdYlBu',#pic
#                    np.arange(0,1.45,0.05))
# plt.colorbar(ax4)
plt.xlabel('$a_0$')
plt.ylabel('$n_0$')

plt.contour(x, y, delta,
                levels = [0.998],
                colors = 'navy')
plt.text(18.6, 6.2, '$W_{pic}/W_{an}=1$', withdash=False,fontdict={'family':'verdana'},size="x-large")
plt.text(13.9, 8,'$p_{cr}/a_0=1$', withdash=False,fontdict={'family':'verdana'},size="x-large")

# plt.savefig('pcr_a03.pdf')


# In[342]:


A_pad = [10,15,20,25]
P_cr = [] 
# n_list = [2,5,8,11,14,17]
# a_list = [10,15,20,25]
n_list = [4,6,8,10]

P = np.zeros((len(A_pad),len(n_list)))
P_A_ = np.zeros((len(A_pad),len(n_list)))


for iii in range(len(n_list)):
    n_0 = n_list[iii] 

    print(a_pad,'a_pad задано')
    for j in range(len(A_pad)): 
        
        m = 1.0/(1836) #\mu
        a_pad = A_pad[j]

        print(n_0,'n_0 задано')
       
        b = a_pad/(sqrt(n_0/m)+a_pad)
        ge = 1/sqrt(1-b**2)
        gi = ge
        print(ge,'gamma')

        p = n_0*b

        #начальное значение
        a_0=0.001

        #граничные точки
        z_0=0
        z_f=2.5

        l=p/sqrt(ge**2-1)

        phi_0=0
        a1_0=sqrt(l-1)*a_0
        phi1_0=2*sqrt(l-1)*phi_0


        yy = np.ndarray(shape=(4,1))

        def f(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        p*((ge+y1)/np.sqrt(((y1+ge)**2-(y3)**2)-1)-(gi-y1*m)/np.sqrt((gi-y1*m)**2-(m*y3)**2-1)),
                        y4,
                        -y3+p*y3/np.sqrt(((y1+ge)**2-(y3)**2)-1) ]

        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        k= 0
        try:
            k=solve_ivp(f,(z_0,z_f),y_0)#, max_step=1e-4)
        except:
            print("error")

        phi_p = k.y[0]
        E_p=k.y[1]
        a_p = k.y[2]
        z_p=k.t



        #ионная область
        z_0i=k.t[-1]
        z_fi=z_f+0.5
        print(z_0i,'граница ионы-плазма')

        phi_0=k.y[0][-1]
        phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
        a_0=k.y[2][-1]
        a1_0=sqrt(4*a_pad**2-a_0**2)
        a1_0_1=k.y[3][-1]
        
        print(a_0,'а на гр ионы-пл')
        print(a_pad,'задан a_pad')

        def f_i(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                    y4,
                    -y3 ]


        y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
        except:
            print("error")


        phi_i = ki.y[0]
        E_i=ki.y[1]
        # print(E_i)
        a_i = ki.y[2]
        z_i=ki.t


        #вакуумная область
        z_0v=ki.t[-1]
        z_fv=z_fi+1.5

        phi_0=ki.y[0][-1]
        phi1_0=ki.y[1][-1]
        a_0=ki.y[2][-1]
        print(a_0,'а_0 на гр вак-ион  1')
        a1_0=ki.y[3][-1]

        def f_v(z,y):
            y1, y2, y3, y4 = y
            return [ y2,
                    0,
                    y4,
                    -y3 ]


        y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

        try:
            kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
        except:
            print("error")

        phi_v = kv.y[0]
        E_v=kv.y[1]
        a_v = kv.y[2]
        z_v=kv.t
        
        a=np.hstack((a_p,a_i,a_v))
        z=np.hstack((z_p,z_i,z_v))
        #plt.plot(z,abs(a))


        E=np.hstack((E_p,E_i,E_v))
        phi=np.hstack((phi_p,phi_i,phi_v))
        #plt.plot(z,E)
        #plt.plot(z,phi)
        
        
        if E_i[-1]<=0:
            print('ddd')
            a=0.1 #точность
            for i in range(0,len(E_i)):
                if E_i[i]>=0:
                    d=i
                # print(d,z[d],E[d])
                    if E_i[i]<E_i[d]:
                        d=i
       

        #print(d,z[d],E[d])
        
        #ионная область
            z_0i=k.t[-1]
            z_fi=z_i[d]
        # print(z_0i,'граница ионы-плазма')

            phi_0=k.y[0][-1]
            phi1_0 = np.sqrt(2*n_0*b/m*np.sqrt((m*phi_0 - ge)**2 - 1))
            a_0=k.y[2][-1]
            a1_0=sqrt(4*a_pad**2-a_0**2)
            a1_0_1=k.y[3][-1]

        #print(a1_0,"а' из связи а и а'")
        #print(a1_0_1, "a' из построенного решения")


            def f_i(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        n_0*b*(-(gi-y1*m)/np.sqrt((gi-y1*m)**2-1)),
                        y4,
                        -y3 ]


            y_0_i = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                ki=solve_ivp(f_i,(z_0i,z_fi),y_0_i, max_step = 1e-4)
            except:
                print("error")


            phi_i = ki.y[0]
            E_i=ki.y[1]
            a_i = ki.y[2]
            z_i=ki.t


        #вакуумная область
            z_0v=ki.t[-1]
            z_fv=z_fi+1.5

            phi_0=ki.y[0][-1]
            phi1_0=ki.y[1][-1]
            a_0=ki.y[2][-1]
        #print(a_0)
            print(a_0,'а_0 на гр вак-ион  2')
        #a_pad=a_0
            a1_0=ki.y[3][-1]

            def f_v(z,y):
                y1, y2, y3, y4 = y
                return [ y2,
                        0,
                        y4,
                        -y3 ]


            y_0 = [phi_0,phi1_0,a_0,a1_0] # начальный вектор

            try:
                kv=solve_ivp(f_v,(z_0v,z_fv),y_0, max_step = 1e-4)
            except:
                print("error")

            phi_v = kv.y[0]
            E_v=kv.y[1]
            a_v = kv.y[2]
            z_v=kv.t
        
            a=np.hstack((a_p,a_i,a_v))
            z=np.hstack((z_p,z_i,z_v))
            #plt.plot(z,abs(a))


            E=np.hstack((E_p,E_i,E_v))
            phi=np.hstack((phi_p,phi_i,phi_v))

        
        # ищем седло
        a=0.005 #точность
        z=np.hstack((z_i,z_v))
        for i in range(0,len(z_i)):
            if ((a_i[i]*ki.y[3][i])/np.sqrt(1+(a_i[i])**2))-E_i[i]<a and ((a_i[i]*ki.y[3][i])/np.sqrt(1+(a_i[i])**2))-E_i[i]>-a:
                print(i)
            
        for i in range(0,len(z_v)):
            if ((a_v[i]*kv.y[3][i])/np.sqrt(1+(a_v[i])**2))-E_v[i]<a and ((a_v[i]*kv.y[3][i])/np.sqrt(1+(a_v[i])**2))-E_v[i]>-a:
        #         print(i,"в вак")
                g=i
        #print(g)  
        
        s=g #индекс, соответствующий седловой точке
        h_cr=np.sqrt(1+(a_v[s])**2)-phi_v[s]
        p_cr=np.sqrt((h_cr+phi_i[0])**2-(a_i[0])**2-1)
        #print(p_cr,'критический импульс электрона на границе с ионами')
        
        P_cr.append(p_cr)
        P[iii][j] = p_cr
        P_A_[iii][j] = p_cr/a_pad


# In[9]:


X,Y = np.meshgrid(A_pad,n_list)


ax = plt.contourf(X,Y,P_A,levels = 100)
#                   np.arange(0,3,0.1))

# plt.colorbar(ax)
plt.contour(X, Y, P_A_,
                levels = [0.8],
                colors = 'green')

# plt.title('$p_{cr}/a_0$')

a_0 = [10,15,20,25]
n_0 = [4,6,8,10]
# plt.legend(('1,','2'))


x,y = np.meshgrid(n_0,a_0)
# ax4 = plt.contourf(x,y,delta,levels = 100)
#                    cmap='spring',
#                    cmap='autumn',
#                    cmap='RdYlBu',#pic
#                    np.arange(0,1.45,0.05))
# plt.colorbar(ax4)
plt.xlabel('$a_0$')
plt.ylabel('$n_0$')

plt.contour(X, Y, delta,
                levels = [0.998],
                colors = 'navy')
# plt.legend(('$1$','$2$'))
plt.text(18.6, 6.2, '$W_{pic}/W_{an}=1$', withdash=False,fontdict={'family':'verdana'},size="x-large")
plt.text(13.9, 8,'$p_{cr}/a_0=1$', withdash=False,fontdict={'family':'verdana'},size="x-large")


# plt.savefig('two_lines08.pdf')

