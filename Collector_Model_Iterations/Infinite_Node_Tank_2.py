# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 13:07:24 2020

@author: maxaj
"""

import numpy as np
import matplotlib.pyplot as PL
import CoolProp.CoolProp as CP
from openpyxl import *

Nodes = 10 #Input the number of nodes to model the tank with

Tfo = 38
Tlo = 35
Ta = 20
Ts_list = [40,39,38,37,36,35,34,33,32,31]
Ts_array = np.array(Ts_list)
Fci_list = [] #Create an empty list for the collector control values to be appended to
Fli_list = [] #Create an empty list for the load control values to be appended to
m_doti_list = [] #Create an empty list for the net flow values to be appended to
Cpi_list = [] #Create an empty list for the node specific heat values to be appended to
Cpfo_list = [] #Create an empty list for collector fluid outlet specific heat values to be appended to
Cplo_list = [] #Create an empty list for load fluid outlet specific heat values to be appended to
delta_T_list = [] #Create an empty list to append the temperature change values to be appended to
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
for ii in range(Nodes):
    
    #Define the fluid entry from the collector control function  
    #Define the statement for the first node
    if ii == 0 and Tfo > Ts_list[ii]:
        Fci = 1
    elif ii == 0 and Tfo < Ts_list[ii]:
        Fci = 0
    #Define the statement for the last node
    elif ii == (Nodes-1) and Tlo < Ts_list[ii]:
        Fci = 1
    #Define the statement for the ith node
    elif Ts_list[ii-1] >= Tfo and Tfo > Ts_list[ii]:
        Fci = 1
    else:
        Fci = 0
    
    #Append the node result to the collector control list
    Fci_list.append(Fci)        
    
    #Define the fluid entry from the load control function
    #Define the statement for the first node
    if ii == 0 and Tlo > Ts_list[ii]:
        Fli = 1
    elif ii == 0 and Tlo < Ts_list[ii]:
        Fli = 0
    #Define the statement for the last node
    elif ii == (Nodes-1) and Tlo < Ts_list[ii]:      
        Fli = 1
    #Define the statement fpr the ith node
    elif Ts_list[ii-1] >= Tlo and Tlo > Ts_list[ii]:
        Fli = 1
    else:
        Fli = 0
        
    Fli_list.append(Fli)

print(Fci_list)
print(Fli_list) 
       
#Write a loop to get the flow rate from node i-1 to node i for every node
for kk in range(Nodes):
    
    #Define a variable to be used to sum all the collector flow contributions
    Fc_tot = 0 
    
    #if kk == 0:
        #Fc_tot = 0
    #else:
    #Write a loop to sum all the collector flow contributions
    for j in range(kk):
        Fc_tot += Fci_list[j]
        
    #print('Fc_tot:',Fc_tot)
    
    #Define a variable to be used to sum all the load flow contributions    
    Fl_tot = 0
    
    #print(kk+1)
    #Write a loop to sum all the load flow contributions
    for j in range(kk,Nodes):
        Fl_tot += Fli_list[j]
        #print(j)
    
    #print('Fl_tot:',Fl_tot)
    
    #Calculate the net flow into node i from node i-1
    #Define the value for node 1 as zero
    if kk == 0:
        m_doti = 0
    #Calculate the flow for the ith node
    else:
        m_doti = m_dotc*Fc_tot - m_dotl*Fl_tot
    
    #Append the flow value to the node flow list  
    m_doti_list.append(m_doti)
    
    #Calculate the specific heat value for the ith node
    Cpi = CP.PropsSI('Cp0mass','T',Ts_list[kk]+273,'P',P,'water')
    
    #Calculate the average temperature between the flow from the collector outlet and the node fluid
    Tfo_av = (Tfo+Ts_list[kk])/2
    #Use the average temperature to calculate an average specific heat coefficient
    Cpfo = CP.PropsSI('Cp0mass','T',Tfo_av+273,'P',P,'water')
    
    #Calculate the average temperature between the flow from the load outlet and the node fluid
    Tlo_av = (Tlo+Ts_list[kk])/2
    #Use the average temperature to calculate an average specific heat coefficient
    Cplo = CP.PropsSI('Cp0mass','T',Tfo_av+273,'P',P,'water')
    
    #Append the node specific heat values to the appropriate lists
    Cpi_list.append(Cpi)
    Cpfo_list.append(Cpfo)
    Cplo_list.append(Cplo)
    
    
    
#print(kk)
print(m_doti_list)
#print(Cpfo_list)
#print(Cplo_list)

#Cpfo = CP.PropsSI('Cp0mass','T',Tfo,'P',P,'water')

#Cplo = CP.PropsSI('Cp0mass','T',Tfo,'P',P,'water')

#Write a loop to calculate the temperature change for every node
for nn in range(Nodes):
    
    #print(nn)
    
    #Calculate the temperature change for the first node
    if nn == 0 and m_doti_list[nn+1] <= 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci_list[nn]*m_dotc*Cpfo_list[nn]*(Tfo-Ts_list[nn]) + Fli_list[nn]*m_dotl*Cplo_list[nn]*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta) + m_doti_list[nn+1]*Cpi_list[nn+1]*(Ts_list[nn]-Ts_list[nn+1]))*time_step
    elif nn == 0 and m_doti_list[nn+1] > 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci_list[nn]*m_dotc*Cpfo_list[nn]*(Tfo-Ts_list[nn]) + Fli_list[nn]*m_dotl*Cplo_list[nn]*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta))*time_step
    #Calculate the temperature change for the final node
    elif nn == (Nodes-1) and m_doti_list[nn] >= 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci_list[nn]*m_dotc*Cpfo_list[nn]*(Tfo-Ts_list[nn]) + Fli_list[nn]*m_dotl*Cplo_list[nn]*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta) + m_doti_list[nn]*Cpi_list[nn-1]*(Ts_list[nn-1]-Ts_list[nn]))*time_step
    elif nn == (Nodes-1) and m_doti_list[nn+1] < 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci_list[nn]*m_dotc*Cpfo_list[nn]*(Tfo-Ts_list[nn]) + Fli_list[nn]*m_dotl*Cplo_list[nn]*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta))*time_step
    #Calculate the temperature change for the ith node
    elif m_doti_list[nn] >= 0 and m_doti_list[nn+1] <= 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci_list[nn]*m_dotc*Cpfo_list[nn]*(Tfo-Ts_list[nn]) + Fli_list[nn]*m_dotl*Cplo_list[nn]*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta) + m_doti_list[nn]*Cpi_list[nn-1]*(Ts_list[nn-1]-Ts_list[nn]) + m_doti_list[nn+1]*Cpi_list[nn+1]*(Ts_list[nn]-Ts_list[nn+1]))*time_step
    elif m_doti_list[nn] >= 0 and m_doti_list[nn+1] > 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci_list[nn]*m_dotc*Cpfo_list[nn]*(Tfo-Ts_list[nn]) + Fli_list[nn]*m_dotl*Cplo_list[nn]*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta) + m_doti_list[nn]*Cpi_list[nn-1]*(Ts_list[nn-1]-Ts_list[nn]))*time_step
    elif m_doti_list[nn] < 0 and m_doti_list[nn+1] <= 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci_list[nn]*m_dotc*Cpfo_list[nn]*(Tfo-Ts_list[nn]) + Fli_list[nn]*m_dotl*Cplo_list[nn]*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta) + m_doti_list[nn+1]*Cpi_list[nn+1]*(Ts_list[nn]-Ts_list[nn+1]))*time_step
    elif m_doti_list[nn] < 0 and m_doti_list[nn+1] > 0:
        delta_T = (1/(mass*Cpi_list[nn]))*(Fci_list[nn]*m_dotc*Cpfo_list[nn]*(Tfo-Ts_list[nn]) + Fli_list[nn]*m_dotl*Cplo_list[nn]*(Tlo-Ts_list[nn]) - U_tank*A_tank*(Ts_list[nn]-Ta))*time_step
    else:
        print('You fucked it up you prick')
    
    #Append the temperature change values to the node temperature change list    
    delta_T_list.append(delta_T)
    
#Convert the node temperature change list into a numpy array
delta_T_array = np.array(delta_T_list)

#Calculate the new node temperatures
Ts_array = Ts_array + delta_T_array

print(delta_T_array)
print(Ts_array)