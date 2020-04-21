# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 11:23:28 2020

@author: maxaj
"""

import numpy as np

Ib = 850 #Input incident beam radiation in watts per metres squared
Id = 150 #Input incident diffuse radiation in watts per metres squared
theta = 0 #Input angle of incidence of radiation on collector in radians
beta = 0 #Input slope of collectors in degrees
Pg = 0.2 #Input value of ground reflectance
t = 0.926 #Input collector transmitance 
a = 0.95 #Input collector absorbance
Ac = 1.84 #Input collector area in metres squared


def absorbed_radiation(Ib,Id,theta,beta):
    
    #Calculate total radiation incident on a flat surface
    #I = Ib+Id
    #print('I:',I)
        
    #Allow function to be exited early if there is no solar radiation incident on the collector
    #if I == 0 and Ibn == 0:
    #Qu = 0
    #It_col = 0
    #Function_outputs = [Qu,It_col]
    #return Function_outputs
        
        
    #Calculate B (in degrees) for use in the delta calculation
    #B_deg = (n-1)*(360/365)
    
    #Convert B into radians
    #B = B_deg*np.pi/180
    
    #Calculate delta (in degrees) - declination - angle of the sun at solar noon i.e. whent the sun is on the local meridian
    #delta_deg = (180/(np.pi))*(0.006918-0.399912*np.cos(B)+0.070257*np.sin(B)-0.006758*np.cos(2*B)+0.000907*np.sin(2*B)-0.002697*np.cos(3*B)+0.00148*np.sin(3*B))
    
    #Convert delta into radians
    #delta = delta_deg*np.pi/180
    
    #Calculate omega (in degrees) - hour angle
    #omega_deg = h*15-187.5
    #Convert omega to radians
    #omega = omega_deg*np.pi/180
    
    #Calculate theta (in radians) - Angle of incidence on the collector plate
    #theta = np.arccos((np.sin(delta)*np.sin(lat-beta))+(np.cos(delta)*np.cos(omega)*np.cos(lat-beta)))
    
    #Convert theta to degrees for use in angular beam transmition absorption product calculation
    theta_deg = theta*180/np.pi
        
    #Allow function to be exited early if there is no solar radiation incident on the collector
    #if I == 0 and Ibn == 0:
        #Qu = 0
        #It_col = 0
        #Function_outputs = [Qu,It_col,theta_deg]
        #return Function_outputs
        
    I = Ib+Id
    
    #Calculate the beam radiation incident on the tilted collector surface
    Ibt = Ib*(np.cos(theta))

    #Calculate the beam radiation incident on a horizontal surface
    #Ibh = Ibn*((np.sin(delta))*(np.sin(lat))+(np.cos(delta))*(np.cos(omega))*(np.cos(lat)))
    
    #Calculate the diffuse radiation incident on a horizontal surface
    Idh = Id
    
    #Calculate equivalent theta for ground reflected radiation
    theta_tag = 90-(0.5788*beta)+(0.002693*(beta**2))
    #print('theta_tag:',theta_tag)
    
    #Calculate equivalent theta for diffuse radiation
    theta_tad = 59.7-(0.1388*beta)+(0.001497*(beta**2))
    #print('theta_tad:',theta_tad)
    
    #Calculate normal transimition absorption product
    tan = 1.01*t*a
    #print('tan:', tan)
    
    #Calculate the angular transmition absorption product for beam radiation - using IAM specific to TVP collector
    if theta_deg <= 30:
        #If theta is less than 30 then transmission absorption product is just equal to the normal value
        tab = tan
    elif theta_deg > 30 and theta_deg < 90:
        #If theta is between 30 and 90 degrees transmission absortion product is given by the polynomial
        tab = tan*((-1*(10**(-8))*theta_deg**4)+(-2*(10**(-6))*theta_deg**3)+(0.0002*theta_deg**2)-(0.003*theta_deg)+1.0061)
    else:
        tab = 0
        print('No beam radiation incident on the collector')
        
    #print('tab:',tab)
    
    #Calculate the angular transmition absorption product for ground reflected radiation
    if theta_tag <= 30:
        #If theta is less than 30 then transmission absorption product is just equal to the normal value
        tag = tan
    elif theta_tag > 30 and theta_tag < 90:
        #If theta is between 30 and 90 degrees transmission absortion product is given by the polynomial
        tag = tan*((-1*(10**(-8))*theta_tag**4)+(-2*(10**(-6))*theta_tag**3)+(0.0002*theta_tag**2)-(0.003*theta_tag)+1.0061)
    else:
        tag = 0
        
    #print('tag:',tag)
    
    #Calculate the angular transmition absorption product for diffuse radiation
    if theta_tad <= 30:
        #If theta is less than 30 then transmission absorption product is just equal to the normal value
        tad = tan
    elif theta_tad > 30 and theta_tad < 90:
        #If theta is between 30 and 90 degrees transmission absortion product is given by the polynomial
        tad = tan*((-1*(10**(-8))*theta_tad**4)+(-2*(10**(-6))*theta_tad**3)+(0.0002*theta_tad**2)-(0.003*theta_tad)+1.0061)
    else:
        tad = 0
        
    #print('tad:', tad)
    
    #Calculate Rb
    #Rb = ((np.cos(lat-beta))*(np.cos(delta))*(np.cos(omega))+(np.sin(lat-beta))*(np.sin(delta)))/((np.cos(lat))*(np.cos(delta))*(np.cos(omega))+(np.sin(lat))*(np.sin(delta)))
    #Rb = 1
    
    #Calculate the absorbed solar radiation per unit collector area
    #S = Ib*Rb*tab + Id*tad*((1+np.cos(beta))/2) + I*Pg*tag*((1-np.cos(beta))/2)
    #print('S:',S)
    
    #Calculate the absorbed solar radiation per unit collector area
    S = Ac*(Ibt*tab + Idh*tad*((1+np.cos(beta))/2) + I*Pg*tag*((1-np.cos(beta))/2))
    #print('S:',S)
    
    return S

print(absorbed_radiation(Ib,Id,theta,beta))