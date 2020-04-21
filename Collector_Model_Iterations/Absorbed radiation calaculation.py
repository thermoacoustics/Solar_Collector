# -*- coding: utf-8 -*-
"""
Created on Sun Mar  1 21:03:50 2020

@author: maxaj
"""

import numpy as np
import matplotlib.pyplot as PL

def absorbed_radiation(I,Ib,Id,Pg,beta,lat,n,h,t,a):
    
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
    
    #Calculate equivalent theta for diffuse radiation
    theta_tad = 90-(0.5788*beta)+(0.002693*(beta**2))
    
    #Calculate equivalent theta for ground reflected radiation
    theta_tag = 59.7-(0.1388*beta)+(0.001497*(beta**2))
    
    #Calculate normal transimition absorption product
    tan = 1.01*t*a
    
    #Calculate the angular transmition absorption product for beam radiation
    tab = theta_deg*tan
    
    #Calculate the angular transmition absorption product for diffuse radiation
    tad = theta_tad*tan
    
    #Calculate the angular transmition absorption product for ground reflected radiation
    tag = theta_tag*tan
    
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
    
    return Sa

    