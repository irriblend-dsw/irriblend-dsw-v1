#!/usr/bin/env python
# coding: utf-8



__author__ = "Belen Gallego-Elvira"
__email__ = "belen.gallego@upct.es"
__version__ = "v1.0" 


########################## CODE STRUCTURE #########################################

# There are four inter-connected parts with specific functionalities
# Input-data import (integrated in fuction "read_input_and_run_simulation")
# Water blend generator (integrated in fuction "read_input_and_run_simulation" 
#                        and calls function "water_blend")
# Fertigation optimization (calls function "ferti_eco_min")
# Visualization and output (integrated in fuction "read_input_and_run_simulation")


###################### IMPORT SCIENTIFIC LIBRARIES ################################

import pandas as pd

import numpy as np

from scipy.optimize import minimize 

from itertools import permutations, \
                      combinations_with_replacement, \
                      product, \
                      combinations

import matplotlib.pyplot as plt

import plotly.graph_objects as go



######################### function "water_blend" ##################################

def water_blend(waters_xlsx, selected_waterss, water_percents):
    """
    blend of water sources
    
    Parameters
    ----------
    waters_xlsx : xlsx file with water ions, CE and price
        see user_input.xlsx, units mmol/l, dS/m, euro/m3
    selected_waterss : list of ints
        selected water sources indexes
    water_percents : list of floats
        percent of each source, e.g. [0.8, 0.1, 0.1]

    Returns  
    -------
    water_blend_CE : float
        CE in dS/m 
    water_blend_price : float 
        prize in euro/m3
    water_blend_mmol_l : list of floats
        [NO3-, NH4+, H2PO4-, K+,  Ca2+,  Mg2+, SO42-,  Cl-, Na, HCO3-] (mmol/l)

    Notes
    -----
    Notes 

    """
    ## Read input data
    df_waters = pd.read_excel(waters_xlsx, sheet_name= 'waters', skiprows = 0,
                       index_col =0, usecols=range(0,13))

    df_waters.fillna(0, inplace=True)

    selected_waters = df_waters.iloc[selected_waterss]

    water_percents = water_percents
    
    if len(water_percents)==2:
        blend = np.sum(selected_waters.values[:,:]*np.array([[water_percents[0]],[water_percents[1]]]), axis = 0)
    elif len(water_percents)==3:
        blend = np.sum(selected_waters.values[:,:]*np.array([[water_percents[0]],[water_percents[1]], [water_percents[2]]]), axis = 0)
    water_blend_price = blend[-1]
    water_blend_CE = blend[-2]

    # molecular_weights : mw
    mw = np.array([62, 97, 96.1, 35.5, 61, 18, 39.1, 40.1, 24.3, 23])
    water_blend_mmol_l = blend[0:-2]/mw
       
    # Order water_blend_ions(mol/l) as   
    #                                      [NO3-, NH4+, H2PO4-, K+,  Ca2+,  Mg2+, SO42-,  Cl-, Na, HCO3-]
    water_blend_mmol_l = water_blend_mmol_l[[0,     5,    1,     6,   7,     8,     2,     3,   9,  4]]
    

    return water_blend_CE, water_blend_price, water_blend_mmol_l


####################### function "ferti_eco_min" ###############################

def ferti_eco_min(ferts_xlsx, ideal_sol, irri_water_blend, irri_water_blend_CE, max_CE):
    """
    Optimized fertilizer selection with scipy.optimize.minimize
    Given a set of commercial fertilizers (user input):
          choose the cheapest combination that provides a given (user input) fertigation solution
    Constrains: CE < CEmax (user input)
                Exact amounts of NO3-, NH4+, H2PO4-, K+, HCO3-
                Amounts of Ca2+, Mg2+, SO42-, Cl-, Na+, equal or above needs.
                 

    Parameters
    ----------
    ferts_xlsx : xlsx file with available fertilizers
        file with available fertilizers, see user_input.xlsx (mmol/g, euro/kg).
    ideal_sol : list with float numbers
        [NO3-, NH4+, H2PO4-, K+,  Ca2+,  Mg2+, SO42-,  Cl-, Na+, HCO3-] (mmol/l)
    irri_water_blend : list with float numbers
        [NO3-, NH4+, H2PO4-, K+,  Ca2+,  Mg2+, SO42-,  Cl-, Na+, HCO3-] (mmol/l)
    irri_water_blend_CE : float
        water blend CE in dS/m
    max_CE: float
        max CE in dS/m (user input)

    Returns  
    -------
    ferti_cost :  float
        fertilizers cost (euro/m3)
    ferti_list : list of strings
        list with the names of the selected fertilizers
    ferti_list_k : list of floats 
        list with the quantities of the selected fertilizers
        (grs of each fertilizer per liter of water)
    sol_final_fertigation : list
        Nutrient solution (water+fertilizers)
        [NO3-, NH4+, H2PO4-, K+,  Ca2+,  Mg2+, SO42-,  Cl-, Na, HCO3-] (mmol/l)


    Notes
    -----
    Notes 

    """
    
    ideal_sol = np.array(ideal_sol)
    irri_water_blend = np.array(irri_water_blend)
    
    
    # Read fertilizer data from user_input.xlsx to a pandas data_frame: df_fert
    df_fert = pd.read_excel(ferts_xlsx, sheet_name= 'ferts', skiprows = 0,
                       index_col =0, usecols=range(0,12))

    df_fert.fillna(0, inplace=True)
    
    
    # Set fertilizer needs
    
    fert_needs = ideal_sol - irri_water_blend
    fert_needs_H = fert_needs[-1]
    fert_needs[fert_needs < 0] = 0
    # fert_needs, negative sign of HC03,  given as H+
    fert_needs[-1] = fert_needs_H
    
    
    # CE constrains
    
    ## CEmax_fert: maximum contribution of ferts to CE
    CEmax_fert = max_CE - irri_water_blend_CE

    ## Anions (meq/l)
    matrix_anions = df_fert.values[:, 0:-1]*np.array([1,0,1,0,0,0,2,1,0,0])

    ## column anions (meq/l) added by each fertilizer
    ## The sum of anions in meq/l divided by 10 is an estimation of CE
    ## ref: American Water Works Association, Water Environment Federation(1999)
    col_anions = np.sum(matrix_anions, axis = 1)
    
    
    # Inputs for scipy.optimize.minimize ('SLSQP' method)
    
    ## Objective: minimize sum(k1*p1+....+kn*pn)
    ## ki: fertiliz amount 
    ## pi: prize 
    
    ## p: objective to be minimized: prize
    p = df_fert.iloc[:,-1].values
    
    ## matrix for eqs and ineqs
    number_fertilizs = np.shape(df_fert.iloc[:,0].values)[0]
    number_elements = 10
    
    ## A_eq: column per element. i.e. amount of element i in each fertilizer
    A_eq = np.transpose(df_fert.values[:, 0:number_elements])

        
    ## CE constrain
    A_ub = np.array([col_anions])/10.
    
    ## b_ub, right side of ineqs
    b_ub  = CEmax_fert
    
   

    def objective(x):
        return np.sum(p*x)


    # NO3-, exact amount in needs  
    def constraint0(x):
        i=0
        return np.sum(x*A_eq[i])-fert_needs[i]
        

    # NH4+, exact amount in needs       
    def constraint1(x):
        i=1
        return np.sum(x*A_eq[i])-fert_needs[i]


    # H2PO4-, exact amount in needs  
    def constraint2(x):
        i=2
        return np.sum(x*A_eq[i])-fert_needs[i]

    # K+, exact amount in needs  
    def constraint3(x):
        i=3
        return np.sum(x*A_eq[i])-fert_needs[i]
    

    # Ca++, must be above need, no max
    def constraint4(x):
        i=4
        return np.sum(x*A_eq[i])-fert_needs[i]
 

    # Mg+,  must be above need, no max
    def constraint5(x):
        i=5
        return np.sum(x*A_eq[i])-fert_needs[i]


    # SO4-2, must be above need, no max
    def constraint6(x):
        i=6
        return np.sum(x*A_eq[i])-fert_needs[i]
        

    # Cl-, must be above need, no max
    def constraint7(x):
        i=7
        return np.sum(x*A_eq[i])-fert_needs[i]


    # Na+, must be above need, no max
    def constraint8(x):
        i=8
        return np.sum(x*A_eq[i])-fert_needs[i]
        
    # H+, exact amount in needs  
    def constraint9(x):
        i=9
        return np.sum(x*A_eq[i])-fert_needs[i]
    
    # CE < CEmax_fert
    # Warning, this may mean that nutrition needs are not satitisfied
    def constraint10(x):
        return CEmax_fert - np.sum(x*np.array([col_anions])/10.)
        


    cons = [{'type':'eq', 'fun': constraint0},
            {'type':'eq', 'fun': constraint1},
            {'type':'eq', 'fun': constraint2},
            {'type':'eq', 'fun': constraint3},
            {'type':'ineq', 'fun': constraint4},
            {'type':'ineq', 'fun': constraint5},
            {'type':'ineq', 'fun': constraint6},
            {'type':'ineq', 'fun': constraint7},
            {'type':'ineq', 'fun': constraint8},
            {'type':'eq', 'fun': constraint9},
            {'type':'ineq', 'fun': constraint10}]
    

    # initial guesses
    x0 = np.zeros(p.size)

    # bounds
    b = (0,3)
    bnds = p.size*((b,))        

    res = minimize(objective,x0,method='SLSQP', bounds=bnds,constraints=cons)

    ## Return
    
    # Prize of optimal fertilizers mix 
    ferti_cost = round(res.fun,4)
    
    # list of selected fertilizers and amounts 
    # (k: grs of each fertiliz to add per liter of water)
    list_fertilizers = list(df_fert.index.values)
    ferti_list = []
    ferti_list_k = []
    for q,f in zip(res.x, list_fertilizers):
        if q>0.00001: # remove traces
            ferti_list.append(f)
            ferti_list_k.append(round(q,4))
            
    
    ## Final solution (water+fertliz) (mmol/l)
    # water ions [NO3-, H2PO4-, K+, Ca2+,  Mg2+, SO42-,  Cl-, HCO3-]
    
    mx_data = df_fert.values[:, 0:10]
    mx_q_fert = mx_data*res.x[:, None]
    sol_fer_calc = np.sum(mx_q_fert, axis=0)    
    sol_final_fertigation = (np.around((irri_water_blend+np.array(sol_fer_calc)),1)).tolist()
    

    return ferti_cost, ferti_list, ferti_list_k, sol_final_fertigation



########### MAIN FUNCTION: READ INPUT, RUN THE SIMULATION, SAVE OUTPUT #############
################# function "read_input_and_run_simulation" #########################

def read_input_and_run_simulation(user_input_xlsx):
    """
    Read input file and run simulation

    Parameters
    ----------
    user_input_xlsx : xlsx file with user input described in USER_GUIDE
    """
    
    ############################ Input-data import #################################
    
    #----------1) Crop ideal solution ---------------------------------------------#

    df = pd.read_excel(user_input_xlsx, sheet_name = 'fert-solution') 

    # Ideal solution (mmol/l) 
    # [NO3-, NH4+, H2PO4-, K+,  Ca2+,  Mg2+, SO42-,  Cl-, Na, HCO3-]

    ideal_sol = df.iloc[0].tolist()[1:-1]

    # EC limit in fertigation solution

    max_CE = df.iloc[0].tolist()[-1]



    #----------2) Available waters ------------------------------------------------#
    # xlsx file with available fertilizers 
    # (ions: mmol/g,CE: dS/m, prices: euro/m3)

    waters_xlsx = user_input_xlsx

    df = pd.read_excel(user_input_xlsx, sheet_name = 'waters') 

    selected_water_names = df['Water-name'].tolist()

    max_percent_waters = (df['Max (%)']/100.).tolist()

    # if there are more than 3 water sources: keep only the first three
    if len(selected_water_names) >3:

        selected_water_names = selected_water_names[0:3]



    #----------3) Available fertilizers -------------------------------------------#
    # xlsx file with available fertilizers 
    # (ions: mmol/g, prices: euro/kg, density:g/cm3)

    ferts_xlsx = user_input_xlsx

    df = pd.read_excel(user_input_xlsx, sheet_name = 'ferts', skiprows = 0)

    available_ferts = df['Unnamed: 0'].tolist()



    #----------4) Profitability indicators ----------------------------------------#

    df = pd.read_excel(user_input_xlsx, sheet_name = 'crop-profit', skiprows = 1) 

    crop_prof_vals = df['input-values'].tolist()

    # Water productivity (kg yield / m3 water used for fertigation)

    wp = crop_prof_vals[0]

    # Land productivity (kg yield / m2 of crop surface)

    kgcrop__m2 = crop_prof_vals[1]

    # Crop market price (euro/kg) 

    euro__kgCropSold = crop_prof_vals[2]

    # Crop salinity threshold (dS / m) above which there is yield loss 

    CE_threshold = crop_prof_vals[3]

    # percent yield loss per dS/m above crop salinity threshold

    yield_decrease_factor = crop_prof_vals[4]

    ## Note: PPI calculations are made per hectare
    ### Total crop surface: 10000 m2

    cropaream2 = 10000
    

    ########################### Water blend generator ##############################
    
    # selected_waters for the water blend
    ## the program rus with 2 o 3 water sources
    
    
    if len(selected_water_names) == 3 : 
    
      selected_waters = [0,1,2]
    
      ## Possible combinations of selected water sources
    
      ### Percentage step (e.g. 0.05, mix water changing proportions in 5% steps)
    
      perc_step = 0.05
    
      ## all permutations with repetition & 
      ### conds: sum ==1 & indicated max percent of each water source
    
      unique_perc = np.arange(0+perc_step, 1, perc_step)
      permuts = product(unique_perc,repeat=len(selected_waters))
      water_percentss = []
      for urp in permuts:
          p = [round(i,2) for i in urp]
          if (round(np.sum(p),2) == 1        and 
              p[0] <= max_percent_waters[0]  and 
              p[1] <= max_percent_waters[1]  and
              p[2] <= max_percent_waters[2]):    
              water_percentss.append(p)
    
    else:
      
      selected_waters = [0,1]
    
      ### Percentage step (e.g. 0.05, mix water changing proportions in 5% steps)
    
      perc_step = 0.05
    
      ## all permutations with repetition & 
      ### conds: sum ==1 & indicated max percent of each water source
    
      unique_perc = np.arange(0+perc_step, 1, perc_step)
      permuts = product(unique_perc,repeat=len(selected_waters))
      water_percentss = []
      for urp in permuts:
          p = [round(i,2) for i in urp]
          if (round(np.sum(p),2) == 1        and 
              p[0] <= max_percent_waters[0]  and 
              p[1] <= max_percent_waters[1]):
              water_percentss.append(p)
    
    
    
    ## number of possible combinations = len(water_percentss)  
    water_combinations = len(water_percentss) 
    
    
    # lists to store simulation results
    ## each element of the list correspond to a unique water blend 
    ## e.g. [0.05, 0.25, 0.70]
    ## 5% watersource1, 25% watersource2, 70% watersource3
    
    water_blend_CEs = []
    water_blend_prices = []
    water_blend_mmol_ls = []
    ferti_costs = []
    ferti_lists = []
    ferti_list_ks = []
    sol_final_fertigations = []
    sol_final_fertigation_CEs = []
    kg_m3_fertliz_mixs = []
    euro_m3_fertliz_mixs = []
    price_m3_solutions = []
    water_percents_with_sols = []
    
    
    # for each combination of water sources generated
    ## 1. Blend water
    ## 2. Find optimize fertigation
    
    for e, water_percents in enumerate(water_percentss):
        # print(e)
        try:
            
            ###### Call water_blend function 
            
            # run water_blend function
            ## water_blend function outputs
            ### water_blend_CE : electrical conductivity of the water blend in dS/m
            ### water_blend_price : cost of water blend in euros/m3 
            ### water_blend_mmol_l : water blend ions (mmol/l) 
            ### [NO3-, NH4+, H2PO4-, K+,  Ca2+,  Mg2+, SO42-,  Cl-, Na+, HCO3-]        
            
            water_blend_CE, water_blend_price, water_blend_mmol_l = water_blend(waters_xlsx, 
                                                                             selected_waters, 
                                                                             water_percents)
     
    
            ########################## Fertigation optimization ###############################
            # run ferti_eco function
            '''
            Given a set of commercial fertilizers (user input), ferti_eco chooses the 
            combination that provides a given (user input) fertigation solution at the lowest cost.
            Constrains: CE < CEmax (user input)
                        Exact amounts of NO3-, NH4+, H2PO4-, K+, HCO3-
                        Amounts of Ca2+, Mg2+, SO42-, Cl-, Na+, equal or above needs
            '''
            ## ferti_eco function outputs
            ### ferti_cost : fertilizers cost (euro/m3)
            ### ferti_list : list with the names of selected fertilizers
            ### ferti_list_k : list with the quantities of selected fertilizers
                        #### (grams of each fertilizer to add per liter of water)
            ### sol_final_fertigation : nutrient solution (water+fertilizs) (mmol/l)
                #### [NO3-, NH4+, H2PO4-, K+,  Ca2+,  Mg2+, SO42-,  Cl-, Na, HCO3-]
            
            ferti_cost, ferti_list, ferti_list_k, sol_final_fertigation = ferti_eco_min(ferts_xlsx, 
                                                                                       ideal_sol, 
                                                                                       water_blend_mmol_l, 
                                                                                       water_blend_CE, 
                                                                                       max_CE) 
    
    
            # estimate CE of fertigation solution as sum(anions)/10 
            ## 100 × anion (or cation) sum, meq/L = (0.9–1.1) EC
            ## Reference: American Public Health Association, 
            ### American Water Works Association, 
            ### Water Environment Federation(1999)
            ### Standard Methods for the Examination of Water and Wastewater
            
            sol_final_fertigation_CE = np.sum(np.array(sol_final_fertigation)*[1,0,1,0,0,0,2,1,0,0])/10.
              
           
            # filter outputs of ferti_eco function without optimal solution      
            if sol_final_fertigation[0] < ideal_sol[0]:
                error_go_to_except
            
            
            # Calculate costs
    
            ## fertilizer cost 
            kg_m3_fertliz_mix = np.sum(np.array(ferti_list_k))
            euro_m3_fertliz_mix = 1*ferti_cost
    
            ## fertigation solution cost (fertilizer+water)
            price_m3_solution = water_blend_price + euro_m3_fertliz_mix 
    
            
            # save results for a given combination of water sources      
            
            water_blend_CEs.append(water_blend_CE)
            water_blend_prices.append(water_blend_price)
            water_blend_mmol_ls.append(water_blend_mmol_l)
            ferti_costs.append(ferti_cost)
            ferti_lists.append(ferti_list)
            ferti_list_ks.append(ferti_list_k)
            sol_final_fertigations.append(sol_final_fertigation)
            sol_final_fertigation_CEs.append(sol_final_fertigation_CE)
            kg_m3_fertliz_mixs.append(kg_m3_fertliz_mix)
            euro_m3_fertliz_mixs.append(euro_m3_fertliz_mix)
            price_m3_solutions.append(price_m3_solution) 
            water_percents_with_sols.append(water_percents)
            
        except:
            excpt = 'optimum fertilizer combination not found for this water blend'
    
    
    #!   The program continues as long as there is at least 1 possible solution   !#
    
    if len(water_percents_with_sols) > 0:
    
      # Yield loss due to salinity 
      ## if CE is over CE_threshold, how is the yield affected?
      ## FAO 61, Annex 1. Crop salt tolerance data: 
      ## yield is reduced linearly over CE_threshold
      ## Relative yield: Yr = 100 - b(ECr - a)
      ## b = the slope expressed in percent per dS/m (= yield_decrease_factor)
      ## a = the salinity threshold expressed in dS/m (= CE_threshold)
      ## ECr = EC in the root zone of the soilless culture in dS/m
    
      reduction_salt_estress = []
      for ce in sol_final_fertigation_CEs :
        
        if ce < CE_threshold:
            reduction_salt_estress.append(1)
            
        else:
            relative_yield = (100-yield_decrease_factor*(ce - CE_threshold))/100        
            reduction_salt_estress.append(relative_yield)
            
              
    
    
      # Potential profitability calculations 
      ## Known variables:
      ### kg/m3  = reduction_salt_estress*water_productivity
      kg_crop__m3_sol = np.array(reduction_salt_estress)*wp
    
      ### euro/m3 = price_m3_solutions
      arr_price_m3_solutions =  np.array(price_m3_solutions)
    
      ### how much does it cost to produce 1 kg crop for each water combination?
      euro__kg_crop = arr_price_m3_solutions/kg_crop__m3_sol
    
      ### index (location) of minimun price per kg
      min_price_loc = np.argmin(euro__kg_crop)
    
      euro__kgCropProduced = euro__kg_crop
    
      # Potential Benefit
    
      cash_in = kgcrop__m2 * cropaream2 * euro__kgCropSold * np.array(reduction_salt_estress)
    
      cash_out = kgcrop__m2 * cropaream2 * euro__kgCropProduced * np.array(reduction_salt_estress)
    
      benefit = cash_in - cash_out
    
      ## combination of waters that provides the maximum benefit
      max_benefit_loc = np.argmax(benefit)
    
    
      ######################## Visualization and output #####################################
    
      # Create a Pandas Excel writer to save output
    
      writer = pd.ExcelWriter('/content/drive/My Drive/irriblend-dsw-data/output.xlsx')
    
      # water combinations
    
      col_names = [selected_water_names[i] for i in selected_waters] + \
                  ['CE (dS/m)']                                      + \
                  ['price (\u20ac/m3)'] 
    
      data_dff = [list(water_percents_with_sols[i]) + 
                [water_blend_CEs[i]]               +
                [water_blend_prices[i]]           
                for i in range(0,len(water_percents_with_sols))]
                
    
      dff = pd.DataFrame(data=data_dff,    
                        columns=col_names)
    
      dff.to_excel(writer, sheet_name='water_combs')
    
    
      # optimum fertilizer combination for a given water combination
    
      col_names = available_ferts
    
      list_fert_ks = []
    
      for i in range(0,len(water_percents_with_sols)):
          
          kfa = np.zeros(len(available_ferts))
    
          indexs=[available_ferts.index(x) for x in ferti_lists[i]]
    
          kfa[indexs] = ferti_list_ks[i]
    
          list_fert_ks.append(kfa.tolist())
    
    
      data_dff = list_fert_ks
                
    
      dff = pd.DataFrame(data=data_dff,    
                        columns=col_names)
    
    
      dff.to_excel(writer, sheet_name='opt_fert_comb')
    
    
      # fertigation solutions
    
      col_names = ['NO3-','NH4+','H2PO4-','K+', 'Ca2+','Mg2+','SO42-','Cl-', 
                  'Na+', 'HCO3-', 'CE dS/m', 'price (\u20ac/m3)']
    
    
      data_dff = [list(sol_final_fertigations[i])   + 
                [sol_final_fertigation_CEs[i]]     +
                [price_m3_solutions[i]]           
                for i in range(0,len(water_percents_with_sols))]
                
    
      dff = pd.DataFrame(data=data_dff,    
                        columns=col_names)
    
    
      dff.to_excel(writer, sheet_name='fertig_sol')
    
    
      # potential profit, water productivity, production cost
    
      col_names = ['PPI (\u20ac/ha)']               +\
                  ['water productivity (kg yield/m3)'] +\
                  ['production cost (\u20ac/kg yield)']
    
    
      data_dff= [[benefit[i]]                      +
                [kg_crop__m3_sol[i]]              +
                [euro__kg_crop[i]]
                for i in range(0,len(water_percents_with_sols))]
    
    
      dff = pd.DataFrame(data=data_dff,    
                        columns=col_names)
    
    
      dff.to_excel(writer, sheet_name='profitability')
    
    
      # Optimum water blend
    
      col_names = [selected_water_names[i] for i in selected_waters]
    
      data_dff = [list(water_percents_with_sols[max_benefit_loc])]
                  
    
      dff = pd.DataFrame(data=data_dff,
                        index =  [max_benefit_loc],  
                        columns=col_names)
    
    
      dff.to_excel(writer, sheet_name='optimun_water_blend')
    
    
      # Close the Pandas Excel writer and output the Excel file.
    
      writer.save()
    
    
      ########### print simulation output highlights 
      # Print results
      np.set_printoptions(precision=3)
      np.set_printoptions(suppress=True)
    
      loc = max_benefit_loc
    
      print('')  
      print('************************* Simulation-output highlights **************************')
      print('')
      print('------------------------- SIMULATED WATER COMBINATIONS --------------------------')
      print('')    
      print('The total number of water combinations (5% steps) is: ', water_combinations)
      print('The total number of water combinations with solution (subjected to user constraints*) is: ', len(water_percents_with_sols))
      print('')
      print('*User input constraints:')
      print('')
      print('The percent availability of each water source is:')
      for i in range(len(selected_water_names)):
          print(selected_water_names[i],':', "{0:.0f}".format(max_percent_waters[i]*100), '%')
      print('') 
      print('The EC (dS/m) in the simulated fertigation solution cannot exceed:', "{0:.1f}".format(max_CE)) 
      print('') 
      print('')   
      print('------------------------- OPTIMAL WATER COMBINATION -----------------------------')
      print('')         
      print('The optimal water combination (subjected to user constraints*) is:')
      for i in range(len(selected_water_names)):
          print(selected_water_names[i],':', "{0:.0f}".format(water_percents_with_sols[loc][i]*100), '%')
      print('')
      print('The water_blend CE (dS/m) is: ', "{0:.3f}".format(water_blend_CEs[loc]))
      print(' ')
      print('The water_blend price (\u20ac/m3) is: ',"{0:.3f}".format(water_blend_prices[loc]) )
      print('')
      print('')  
      print('------------- FERTILIZERS FOR THE OPTIMAL WATER COMBINATION ---------------------')   
      print('')    
      print('The fertilizer combination with the lowest cost is the following (name : amount (g/l))')
      for f,k in zip(ferti_lists[loc], ferti_list_ks[loc]):
          print(f,':',k)
    
      print('')    
      print('The prize of the fertilizer (\u20ac/m3) is: ', "{0:.3f}".format(ferti_costs[loc]))
      print('')
      print('')
      print('----------- NUTRIENT SOLUTION FOR THE OPTIMAL WATER COMBINATION -------------------')    
      print('') 
      print('The USER INPUT ideal nutrient solution composition (mmol/l) was:')
      print('[NO3-,NH4+,H2PO4-,K+, Ca2+,Mg2+,SO42-,Cl-, Na, HCO3-]')
      print(ideal_sol) 
      print('')
      print('') 
      print('The nutrient solution composition (mmol/l) computed for the optimal water combination is:')
      print('[NO3-,NH4+,H2PO4-,K+, Ca2+,Mg2+,SO42-,Cl-, Na, HCO3-]')
      print(sol_final_fertigations[loc])
      print('')
      print('The CE (dS/m) of the computed nutrient solution is: ', "{0:.2f}".format(sol_final_fertigation_CEs[loc]))
      print('')
      print('')
      print('--------- COSTS AND PROFIT POTENTIAL INDICATOR (PPI) FOR THE OPTIMAL WATER COMBINATION ----------------')
      print('')  
      print( 'The cost of water+fertilizer(\u20ac/m3) is: ', "{0:.3f}".format(price_m3_solutions[loc]))   
      print('')
      print('The water productivity is (kg_crop/m3):', "{0:.2f}".format(kg_crop__m3_sol[loc]))
      print('')
      print('The fertigation production cost (\u20ac/kg) is: ',  "{0:.3f}".format(euro__kg_crop[loc]))
      print('')
      print('The PPI* (\u20ac/ha) is: ', "{0:.2f}".format(benefit[loc]))
      print('')
      print('*PPI (1ha) =        yield income                   -        fertigation cost')
      print('           (production (kg) x market price (€/kg)) - (cost of water plus fertilizer)')
      print('')
      print('')
      print('')
    
    
      ######### 3D Plot: Profit vs water percent of two waters 
    
      def plot_profit_surface(x_axis, y_axis):
          """
          3D interactive plot of
          Profit (z) vs water percents of waters x and y
          
          Parameters
          ----------
          x_axis : str
              name of water to plot in x axis
          y_axis : str
              name of water to plot in y axis
    
          Returns  
          -------
          None
    
          Notes
          -----
          go.Mesh3d draws a 3D set of triangles with vertices given by x, y and z
          see documentation at https://plotly.com/python/3d-mesh/
          """
    
          # Values for x, y, z axis
    
          x = np.array([num[selected_water_names.index(x_axis)] for num in water_percents_with_sols])*100
    
          y = np.array([num[selected_water_names.index(y_axis)] for num in water_percents_with_sols])*100
    
          z = np.array(benefit)
    
    
          # zz: Max benefit
          # xx, yy: values of x and y at zz
    
          xx = water_percents_with_sols[loc][selected_water_names.index(x_axis)]*100
    
          yy = water_percents_with_sols[loc][selected_water_names.index(y_axis)]*100
    
          zz = benefit[loc]
    
    
          # Plot interactive 3D mesh with Plotly library
          fig = go.Figure(data=[
                        go.Mesh3d(x=x,
                        y=y,
                        z=z,
                        intensity=z,
                        colorscale='viridis',
                        
                        colorbar=dict(title='PPI (\u20ac/ha)', len=0.7),
                        showscale = True),
    
                        go.Scatter3d(name="PPI", 
                                      x=[xx], y=[yy], z=[zz], mode='markers',
                                      marker=dict(size=6))
                                                              
                        ])
    
          fig.update_layout(
              scene = dict(
                      xaxis = dict(nticks=10, ticks='outside',range=[x.min(),x.max()],tickfont=dict(size=15)),
                      yaxis = dict(nticks=10, ticks='outside',range=[y.min(),y.max()],tickfont=dict(size=15)),
                      zaxis = dict(nticks=10, ticks='outside',range=[z.min(),z.max()],tickfont=dict(size=15)),
                      
                      xaxis_title = x_axis+ ' (%)',
                      yaxis_title = y_axis+ ' (%)',
                      zaxis_title='PPI (\u20ac/ha)'
                      
                      ),
              width=1200,
              height=900,
          
              font=dict(size = 15),     
          
              margin=dict(r=20, l=10, b=50, t=10))
              
          fig.update_traces(hovertemplate = x_axis+': %{x}% <br>'+y_axis+': %{y}% <br>' + 'PPI (\u20ac/ha): ' + '%{z:.0f}<extra></extra>')
    
          fig.show()
    
    
      for combi in list(combinations(range(len(selected_waters)),2)):
    
        # Select x axis
        x_axis = selected_water_names[combi[0]]
    
        # Select y axis
        y_axis = selected_water_names[combi[1]]
    
        plot_profit_surface(x_axis, y_axis)
    
    else:
      print('irriblend can not find any solution for the given input data')

    




























