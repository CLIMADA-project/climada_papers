# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 14:00:01 2023

Get damage, and exposure information at hazard centroid for each date. 
The resulting intermediate files are the basis for most figures and stored in a subfolder data/data_at_centroid

how to use:

python get_per_gridcell_data.py -croptypes wheat barley maize rapeseed grapevine -hazard_var MESHS -res 1 2 4 8 16 32

@author: Raphael Portmann
"""
import sys
import os
from climada.engine import Impact
from utility import aggregate_hazard, en_names
from calibration import empirical_calibration_per_exposure, DROP_DATES_DICT
import pickle
import re
import argparse
import warnings 
from climada.util.api_client import Client

def drop_events(drop_dict, damages):
    """
    Drop events from damage impact objects for different croptypes


    Parameters
    ----------
    drop_dict: dict
        Dictionary with events to be dropped per croptype
    damages : dict
        Dictionary of damages (damages as impact objects) from which to drop events

    Returns
    -------
    damages : dict
        Dictionary of damages with dropped events

    """
    for croptype in damages.keys():
        events=list(damages[croptype].event_name)
        for ev in drop_dict[croptype]:
            if ev in events:
                events.remove(ev)
        damages[croptype]=damages[croptype].select(event_names=events)
    return damages

def main(croptypes, hazard_vars, resolutions):  
         
    print('-------------------------------')
    print('GET AT CENTROID DATA:')
    print('-------------------------------')
    print('For')
    print(f'croptypes: {croptypes}')
    print(f'hazard_var: {hazard_vars}')
    print(f'resolutions: {resolutions}')
    print('-------------------------------')

    cwd=os.getcwd()
    outdir=cwd+'/data/data_at_centroid/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f'create {outdir}')

    print(f'Output directory: {outdir}')

    #silence user warnings
    warnings.filterwarnings(action='ignore', category=UserWarning) # setting ignore as a parameter and further adding category    
  
    resolutions=list(resolutions)
    croptypes=list(croptypes)
    hazard_vars=list(hazard_vars)
    
    # read standard data
    print('reading damage and exposure data...')
    client=Client()
    exposure={}
    damages={}
    for croptype in croptypes:
        exposure[croptype] = client.get_exposures('crops', properties={'res_km': '1', 'crop': croptype})
        ds = client.get_dataset_info(data_type='hail_damage_crops', properties={'res_km': '1','crop': croptype})
        _, [impact_file] = client.download_dataset(ds)
 
        damages[croptype] = Impact.from_hdf5(impact_file)
            
    #drop damages at dates where already most of crops were harvested
    damages=drop_events(DROP_DATES_DICT, damages)

    for hazard_var in hazard_vars:
        print(f'GET DATA FOR HAZARD: {hazard_var}')
        #read aggregated hazard data
        hazard_rad=client.get_hazard('hail', properties={'variable': hazard_var})  
        resolutions_agg=resolutions.copy()
        if '1' in resolutions:
            resolutions_agg.remove('1')
        haz={}
        for res in resolutions_agg:
            haz[f'{res}'], _ = aggregate_hazard(hazard_rad, original_grid_epsg = 2056, 
                                  extent_new = None, 
                                  cell_size_new = int(res)*1000,
                                  projection_new_epsg = 4326,
                                  aggfunc = 'max')
        
        #store 1km hazard in haz object
        haz['1']=hazard_rad
    
    for hazard_var in hazard_vars:
        print(f'GET AT CENTROID DATA for: {hazard_var}')
        #create dictionary for output
        at_centroid_data = {}
        for croptype in croptypes:
            print(f'...and croptype: {croptype}')
            
            for j,res in enumerate(resolutions):
                print(f'...and resolution: {res}km')
     
                _, _, _, _, values_at_centroid_all,_ =\
                    empirical_calibration_per_exposure(hazard_object = haz[res],
                                                       exposure_object = exposure[croptype],
                                                       damages = damages[croptype],
                                                       exposure_type = 'agriculture',
                                                       variable = hazard_var,
                                                       fraction_insured = 1)
            
                at_centroid_data[f'{res}km'] = values_at_centroid_all
            
            name = f'data_at_centroid_{hazard_var}_{croptype}.p'
            # save dictionary to pickle file
            with open(outdir+name, 'wb') as file:
                pickle.dump(at_centroid_data, file, protocol=pickle.HIGHEST_PROTOCOL)
                print(f'File {file} saved.')
                    
#setup argument parser and run main
parser = argparse.ArgumentParser()
parser.add_argument(
            "-croptypes", "--croptype",
            dest    = "croptypes",
            nargs   = '+',
            default = ['wheat','barley','maize','rapeseed','grapevine'],
            )
parser.add_argument(
            "-hazard_var", "--hazard variables",
            dest    = "hazard_vars",
            nargs   = "+",
            default = ['MESHS','POH'],
          )     
parser.add_argument(
            "-res", "--resolutions in km(can be one or several of 1, 2, 4, 8, 16, 32",
            dest    = "resolutions",
            nargs   = '+',
            default = ['1','2','4','8','16','33'],
          )       

args = parser.parse_args()   
main(**vars(args))