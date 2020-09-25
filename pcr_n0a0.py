#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np 
from multiprocessing import Process
from math import sqrt
from scipy.integrate import solve_ivp

m_e = 9.1/10**31
e = 1.6/10**19
c = 3*10**10
m_i = m_e*1836


# In[6]:


#задаем n_0 и a_0 в циклах, из них считаем gamma
delta = 4
x_0 = 2
n = 5 # количество точек
A_pad=[25,30,35,40]
P_cr=[] # по z при определенной концентрации
n_list=[15,20,25,30]
P=np.zeros((len(A_pad),len(n_list)))
P_A=np.zeros((len(A_pad),len(n_list)))


for iii in range(len(A_pad)):
    a_pad = A_pad[iii] #концентрация, которая по х
    print(a_pad,'a_0 задано')
    for j in range(len(n_list)):
#         print(j,'j')
        #ЗАДАВАЕМЫЕ ЗНАЧЕНИЯ
        # свободные параметры    
       
        n_0 = n_list[j]
        print(n_0,'n_0 задано')
        m = 1.0/(1836)/5 #\mu
       
    
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
            #plt.plot(z,E)
            #plt.plot(z,phi)
        
        
        
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
        
        #A_pad.append(a_pad)
        P_cr.append(p_cr)
        P[iii][j]=p_cr
        P_A[iii][j]=p_cr/a_pad
#         A_pad.append(a_pad)
print(A_pad)
x,y=np.meshgrid(A_pad,n_list)
print(P_cr)

a=max(P[-1])
b=min(P[0])

d=40

#lev=[]

#  Задаем значение каждого уровня:
#for i in range(10):
 #   lev.append(d)
  #  d=d+3
    

#  Создаем массив RGB цветов каждой области:
#color_region = np.zeros((10, 3))
#color_region[0:, 1:] = 0.2
#color_region[:, 0] = np.linspace(0, 1, 10)


ax1=plt.contourf(x,y,P,cmap='inferno', 
                levels = np.arange(30,70,1))
plt.colorbar(ax1)
plt.xlabel('$a_0$')
plt.ylabel('$n_0$')
#plt.show()

plt.savefig('for_diplom_5.pdf')

