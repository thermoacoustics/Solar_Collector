# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 19:38:21 2020

@author: maxaj
"""

ta = 30
tc = 24
th = 80

  
def table_read(ta,tc,th):
    import pandas as pd
    import numpy as np
    from scipy.interpolate import RegularGridInterpolator
    
    file = r'data_short.xlsx'
    lt = pd.read_excel(file)
    
    c_to_k = 273.15
    
    ta +=c_to_k
    tc +=c_to_k
    th +=c_to_k
    
    ahx_interest = []
    chx_interest = []
    hhx_interest = []
    
    # AHX
    ahx_list     = lt['AHX'] # locate the AHX range and values
    ahx_list     = ahx_list.drop_duplicates(keep='first')
    ahx_interest = ahx_list.iloc[(ahx_list-ta).abs().argsort()[:2]].sort_index().values.tolist() # find the two closest AHX values
    fl_ahx_outrange = (ta > max(ahx_interest) or ta < min(ahx_interest))
    rows = lt[lt['AHX'].isin(ahx_interest)]
        
    
    # CHX
    chx_list     = rows['CHX'] # locate the AHX range and values
    chx_list     = chx_list.drop_duplicates(keep='first')
    chx_interest = chx_list.iloc[(chx_list-tc).abs().argsort()[:2]].sort_index().values.tolist()
    fl_chx_hi    = (tc > max(chx_interest))
    fl_chx_lo    = (tc < min(chx_interest))   
    rows         = rows[rows['CHX'].isin(chx_interest)]
    
    # HHX
    hhx_list     = rows['HHX'] # locate the AHX range and values
    hhx_list     = hhx_list.drop_duplicates(keep='first')
    hhx_interest = hhx_list.iloc[(hhx_list-th).abs().argsort()[:2]].sort_index().values.tolist()
    fl_hhx_hi    = (th > max(hhx_interest))
    fl_hhx_lo    = (th < min(hhx_interest))
    rows         = rows[rows['HHX'].isin(hhx_interest)]
    
    # use the return flags for output
    # when not onset
    if fl_hhx_lo or fl_chx_hi:
        if fl_hhx_lo:
            print('HHX =', th, 'too low, minimum onset', min(hhx_interest))
        elif fl_chx_hi:
            print('CHX =', th, 'too high, maximum onset', max(chx_interest))
        return [0.,0.]
    elif fl_chx_lo or fl_hhx_hi or fl_ahx_outrange:
        print('Undefined temperature range')
        return [0.,0.]
    else:  
        try:
            # transform the 2x2x2 matrix to a interpolatble format qh_mat(ta,tc,th)
            qh_mat = np.zeros((len(ahx_interest),len(chx_interest),len(hhx_interest)))
            qc_mat = np.zeros((len(ahx_interest),len(chx_interest),len(hhx_interest)))
            
            for i in range(len(ahx_interest)):
                for j in range(len(chx_interest)):
                    for k in range(len(hhx_interest)):
                        qh_mat[i,j,k] = rows[(rows['AHX']==ahx_interest[i]) & (rows['CHX']==chx_interest[j]) & (rows['HHX']==hhx_interest[k])]['Q_HHX'].to_numpy()
                        qc_mat[i,j,k] = rows[(rows['AHX']==ahx_interest[i]) & (rows['CHX']==chx_interest[j]) & (rows['HHX']==hhx_interest[k])]['Q_CHX'].to_numpy()
            
            qh = RegularGridInterpolator((ahx_interest,chx_interest,hhx_interest), qh_mat)(np.array([ta,tc,th]))[0]
            qc = RegularGridInterpolator((ahx_interest,chx_interest,hhx_interest), qc_mat)(np.array([ta,tc,th]))[0]
        
            return [qc,qh]
        except ValueError:
            print ('Not onset')
            return [0.,0.]
        
print(table_read(ta,tc,th))