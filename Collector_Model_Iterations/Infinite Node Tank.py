# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:23:57 2020

@author: maxaj
"""

import numpy as np
import matplotlib.pyplot as PL
import CoolProp.CoolProp as CP
from openpyxl import *

Nodes = 10 #Input the number of nodes to model the tank with

Tfo = 50
Tlo = 35
Ta = 20
Ts_list = [40,39,38,37,36,35,34,33,32,31,30]
Ts_array = np.array(Ts_list)
Fci_list = [] #Create an empty list for the collector control values to be appended to
Fli_list = [] #Create an empty list for the load control values to be appended to
m_doti_list = [] #Create an empty list for the net flow values to be appended to
Cpi_list = [] #Create an empty list for the node specific heat values to be appended to
m_dotc = 5 #Input the flow rate through the collector in kg/s
m_dotl = 5 #Input the flow rate through the load in kg/s
P = 10*(10**5) #Input the system pressure in pascals
time_step = 10 #Input the time step for the simulation to be carried out at in seconds

U_tank = 2 #Input the tank loss coefficient in W/m2K
A_tank = 10 #Input the tank surface area in metres squared
V = 2 #Input the tank volume
rho = 1000

mass = V*rho

#Write loop to obtain the collector and load control value for each node
for ii in range(1,Nodes+1):
    
    #Define the fluid entry from the collector control function
    #if ii == 0:
        #Fci = 0    
    if ii == 1 and Tfo > Ts_list[1]:
        Fci = 1
    elif Ts_list[ii-1] >= Tfo and Tfo > Ts_list[ii]:
        Fci = 1
    else:
        Fci = 0
        
    Fci_list.append(Fci)        
    
    #if ii == 0:
        #Fli = 0
    if ii == Nodes and Tlo < Ts_list[Nodes]:
        Fli = 1
    elif Ts_list[ii] >= Tlo and Tlo > Ts_list[ii+1]:
        Fli = 1
    else:
        Fli = 0
        
    Fli_list.append(Fli)

print(Fci_list)
print(Fli_list) 
       
#Write a loop to get the flow rate from node i-1 to node i for every node
for kk in range(1, Nodes+1):
    
    Fc_tot = 0 
    
    if kk == 1:
        Fc_tot = 0
    else:
        for j in range(kk-1):
            Fc_tot += Fci_list[j]
        
    print('Fc_tot:',Fc_tot)
        
    Fl_tot = 0
        
    for j in range(kk,Nodes):
        Fl_tot += Fli_list[j]
    
    print('Fl_tot:',Fl_tot)
    
    if kk == 1:
        m_doti = 0
    else:
        m_doti = m_dotc*Fc_tot - m_dotl*Fl_tot
        
    m_doti_list.append(m_doti)
    
    Cpi = CP.PropsSI('Cp0mass','T',Ts_list[kk]+273,'P',P,'water')
    

print(m_doti_list)

Cpfo = CP.PropsSI('Cp0mass','T',Tfo,'P',P,'water')

Cplo = CP.PropsSI('Cp0mass','T',Tfo,'P',P,'water')

#Write a loop to calculate the temperature change for every node
for nn in range(Nodes):
    
    if m_doti_list[nn] >= 0 and m_doti_list[nn+1] <= 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci[nn]*m_dotc*Cpfo*(Tfo-Ts_list[nn]) + Fli[nn]*m_dotl*Cplo*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta) + m_doti_list[nn]*Cpi_list[nn-1]*(Ts_list[nn-1]-Ts_list[nn]) + m_doti_list[nn+1]*Cpi_list[nn+1]*(Ts_list[nn]-Ts_list[nn+1]))*time_step
    elif m_doti_list[nn] >= 0 and m_doti_list[nn+1] > 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci[nn]*m_dotc*Cpfo*(Tfo-Ts_list[nn]) + Fli[nn]*m_dotl*Cplo*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta) + m_doti_list[nn]*Cpi_list[nn-1]*(Ts_list[nn-1]-Ts_list[nn]))*time_step
    elif m_doti_list[nn] < 0 and m_doti_list[nn+1] <= 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci[nn]*m_dotc*Cpfo*(Tfo-Ts_list[nn]) + Fli[nn]*m_dotl*Cplo*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta) + m_doti_list[nn+1]*Cpi_list[nn+1]*(Ts_list[nn]-Ts_list[nn+1]))*time_step
    elif m_doti_list[nn] < 0 and m_doti_list[nn+1] > 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci[nn]*m_dotc*Cpfo*(Tfo-Ts_list[nn]) + Fli[nn]*m_dotl*Cplo*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta) + m_doti_list[nn+1]*Cpi_list[nn+1]*(Ts_list[nn]-Ts_list[nn+1]))*time_step
    