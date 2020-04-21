# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 17:25:28 2020

@author: maxaj
"""

import numpy as np
import matplotlib.pyplot as PL
import CoolProp.CoolProp as CP

Ib = 42 #Input beam radiation incident on a flat surface in Watts per metre squared
Id = 36 #Input diffuse radiation incident on a flat surface in Watts per metre squared
Pg = 0.5 #Input ground reflected radiation factor
beta = 10 #Input collector slope in degrees
lat = 10 #Input lattitude of the collector location in degrees
n = 1 #Input day of the year (1st January = 1)
h = 12 #Input hour of the day (00:00-01:00 = 1)
t = 0.948 #Input collector transmittance
a = 0.87 #Input collector absorbance
Ac = 1.84 #Input collector area in metres squared
Ul = 5.2 #Input total heat transfer coefficient to the environment in Watts per metres squared kelvin
Tfi = 40 #Input the inlet temperature to the collector in celsius
Ta = 30 #Input the ambient temperature in celsius
m_dot = 0.5 #Input the mass flow rate to the collector in kg per second
F_dash = 0.968 #Input collector efficiency factor
P = 16*(10**5) #Input collector fluid pressure in bar


def absorbed_radiation(Ib,Id,Pg,beta,lat,n,h,t,a,Ac,Ul,Tfi,Ta,m_dot,F_dash,P):
    
    #Calculate total radiation incident on a flat surface
    I = Ib+Id
    
    #Calculate B (in degrees) for use in the delta calculation
    B_deg = (n-1)*(360/365)
    #Convert B into radians
    B = B_deg*np.pi/180
    
    #Calculate delta (in degrees) - declination - angle of the sun at solar noon i.e. whent the sun is on the local meridian
    delta_deg = (180/(np.pi))*(0.006918-0.399912*np.cos(B)+0.070257*np.sin(B)-0.006758*np.cos(2*B)+0.000907*np.sin(2*B)-0.002697*np.cos(3*B)+0.00148*np.sin(3*B))
    #Convert delta into radians
    delta = delta_deg*np.pi/180
    
    #Calculate omega (in degrees) - hour angle
    omega_deg = h*15-187.5
    #Convert omega to radians
    omega = omega_deg*np.pi/180
    
    #Calculate theta (in radians) - Angle of incidence on the collector plate
    theta = np.arccos((np.sin(delta)*np.sin(lat-beta))+(np.cos(delta)*np.cos(omega)*np.cos(lat-beta)))
    
    #Convert theta to degrees for use in angular beam transmition absorption product calculation
    theta_deg = theta*180/np.pi
    
    #Calculate equivalent theta for ground reflected radiation
    theta_tag = 90-(0.5788*beta)+(0.002693*(beta**2))
    
    #Calculate equivalent theta for diffuse radiation
    theta_tad = 59.7-(0.1388*beta)+(0.001497*(beta**2))
    
    #Calculate normal transimition absorption product
    tan = 1.01*t*a
    
    #Calculate the angular transmition absorption product for beam radiation
    if theta_deg <= 30:
        #If theta is less than 30 then transmission absorption product is just equal to the normal value
        tab = tan
    elif theta_deg > 30 and theta_deg < 90:
        #If theta is between 30 and 90 degrees transmission absortion product is given by the polynomial
        tab = tan*((-3*(10**(-6))*theta_deg**3)+(0.0002*theta_deg**2)-(0.0036*theta_deg)+0.9974)
    else:
        tab = 0
        print('No beam radiation incident on the collector')
    
    #Calculate the angular transmition absorption product for ground reflected radiation
    if theta_tag <= 30:
        #If theta is less than 30 then transmission absorption product is just equal to the normal value
        tag = tan
    elif theta_tag > 30 and theta_tag < 90:
        #If theta is between 30 and 90 degrees transmission absortion product is given by the polynomial
        tag = tan*((-3*(10**(-6))*theta_tag**3)+(0.0002*theta_tag**2)-(0.0036*theta_tag)+0.9974)
    else:
        tag = 0
    
    #Calculate the angular transmition absorption product for diffuse radiation
    if theta_tad <= 30:
        #If theta is less than 30 then transmission absorption product is just equal to the normal value
        tad = tan
    elif theta_tad > 30 and theta_tad < 90:
        #If theta is between 30 and 90 degrees transmission absortion product is given by the polynomial
        tad = tan*((-3*(10**(-6))*theta_tad**3)+(0.0002*theta_tad**2)-(0.0036*theta_tad)+0.9974)
    else:
        tad = 0
    
    #Calculate Rb
    Rb = ((np.cos(lat-beta))*(np.cos(delta))*(np.cos(omega))+(np.sin(lat-beta))*(np.sin(delta)))/((np.cos(lat))*(np.cos(delta))*(np.cos(omega))+(np.sin(lat))*(np.sin(delta)))
    
    #Calculate the absorbed solar radiation per unit collector area
    S = Ib*Rb*tab + Id*tad*((1+np.cos(beta))/2) + I*Pg*tag*((1-np.cos(beta))/2)
    
    #Include factor for reduced absorption due to dust
    Df = 0.99
    
    #Include factor for reduced absoprtion due to self shading
    Ds = 0.99
    
    #Calculated adjusted shading when dust and self shading are taken into account
    Sa = S*Df*Ds
    
    #Calaculate initial Cp estimate
    Cp = CP.PropsSI('Cp0mass','T',Tfi+273,'P',P,'water')
    
    #Calculate initial Fr estimate
    Fr = ((m_dot*Cp)/(Ac*Ul))*(1-np.exp(-(Ac*Ul*F_dash)/(m_dot*Cp)))
    
    #Calculate initial estimate of the useful energy gain
    Qu = Ac*Fr*(Sa-Ul*(Tfi-Ta))
    
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
        Cp = CP.PropsSI('Cp0mass','T',Tfm+273,'P',P,'water')
        
        #Calculate new Fr value with new Cp value
        Fr = ((m_dot*Cp)/(Ac*Ul))*(1-np.exp(-(Ac*Ul*F_dash)/(m_dot*Cp)))
        
        #Calculate Tpm estimate
        Tpm = Tfi+(Qu/(Ac*Fr*Ul))*(1-Fr)
        
        #Introduce Qu_est so the value of Qu in this iteration can be compared to the value in the last iteration.
        Qu_est = Qu
        
        #Calculate new Qu estimate
        Qu = Ac*(Sa-Ul*(Tpm-Ta))
        
        #Terminate loop if the new value of Qu is within 1% of the previous value
        if Qu > 0.99*Qu_est and Qu < 1.01*Qu_est:
            break
        else:
            continue
        
    #Calculate the fluid outlet temperature
    Tfo = Ta + (Sa/Ul) + (Tfi-Ta-(Sa/Ul))*np.exp(-(Ac*Ul*F_dash)/(m_dot*Cp))
    
    Qu_check = m_dot*Cp*(Tfo-Tfi)
    
    return Qu, Tfo, Qu_check

print(absorbed_radiation(Ib,Id,Pg,beta,lat,n,h,t,a,Ac,Ul,Tfi,Ta,m_dot,F_dash,P))