# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 15:59:10 2020

@author: maxaj
"""

import numpy as np
import matplotlib.pyplot as PL
import CoolProp.CoolProp as CP

Ac = 1.84 #Input collector area in metres squared
Ul = 5.2 #Input total heat transfer coefficient to the environment in Watts per metres squared kelvin
Tfi = 

def calculate_energy_gain(Ac,Ul,Tfi,Ta,m_dot,F_dash,P,S):
    
    #Calaculate initial Cp estimate
    Cp = CP.PropsSI('Cp0mass','T',Tfi,'P',P,'water')
    
    #Calculate initial Fr estimate
    Fr = ((m_dot*Cp)/(Ac*Ul))*(1-np.exp(-(Ac*Ul*F_dash)/(m_dot*Cp)))
    
    #Calculate initial estimate of the useful energy gain
    Qu = Ac*Fr*(S-Ul*(Tfi-Ta))
    
    #Define i so that a counter for the loop can be used
    i = 2
    Z = 100
    
    for n in range(Z):
        
        #Use counter to inform user on whether or not the iteration has converged to within 1%
        i += 1
        if i == Z:
            print('Qu iteration has not converged to within 1%')
        
        #Calculate new Tfm value
        Tfm = Tfi+(Qu/(Ac*Fr*Ul))*(1-(Fr/F_dash))
        
        #Calculate new Cp value using Tfm
        Cp = CP.PropsSI('Cp','T',Tfm,'P',P,'water')
        
        #Calculate new Fr value with new Cp value
        Fr = ((m_dot*Cp)/(Ac*Ul))*(1-np.exp(-(Ac*Ul*F_dash)/(m_dot*Cp)))
        
        #Calculate Tpm estimate
        Tpm = Tfi+(Qu/(Ac*Fr*Ul))*(1-Fr)
        
        #Introduce Qu_est so the value of Qu in this iteration can be compared to the value in the last iteration.
        Qu_est = Qu
        
        #Calculate new Qu estimate
        Qu = Ac*(S-Ul*(Tpm-Ta))
        
        #Terminate loop if the new value of Qu is within 1% of the previous value
        if Qu > 0.99*Qu_est and Qu < 1.01*Qu_est:
            break
        else:
            continue
        
    #Calculate the fluid outlet temperature
    Tfo = Ta + (S/Ul) + (Tfi-Ta-(S/Ul))*np.exp(-(Ac*Ul*F_dash)/(m_dot*Cp))
    
    return Qu Tfo
        
    
        

        
    