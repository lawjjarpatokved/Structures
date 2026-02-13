import numpy as np
import math as math
import pandas as pd
import openseespy.opensees as ops
import os
import opsvis as opsv
from Structures_2D import Structures_2D
from libdenavit.OpenSees.plotting import *
from libdenavit.OpenSees.get_fiber_data import *
from Units import *
from math import pi, ceil
from libdenavit.section.wide_flange import *
from libdenavit.OpenSees import AnalysisResults
from libdenavit import find_limit_point_in_list, interpolate_list
import copy
import matplotlib.cm as color
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation
import inspect

#################################################
density_of_steel=7850*kg/(m**3)
g=9.81*m/(sec**2)
E=29000*ksi
G=77221*Mpa
Fy=36*ksi
Hk = 0.001*E           # Kinematic hardening modulus

class Steel_Material:
    def __init__(self,mat_tag,E,Fy,G,Hk,density):
        self.mat_tag=mat_tag
        self.E=E
        self.Fy=Fy
        self.G=G
        self.Hk=Hk
        self.density=density
        self.b=Hk / (E + Hk) 

class Moment_Frame_2D(Structures_2D):

    print_ops_status=True

    def __init__(self,width_of_bay,storey_height,
                no_of_elements_column, no_of_elements_beam,
                 beam_section,column_section,load_combination_multipliers,Frame_id,
                 **kwargs):
        # ---- Save a "constructor snapshot" for later cloning ---- This is useful for resetting or duplicating the model. Eg. Calculation of del2_over_del1
        self._init_spec = copy.deepcopy({k: v for k, v in locals().items() if k != "self"})

        self.bay_width=[0]+width_of_bay
        self.storey_height = [0] + storey_height
        self.length_of_frame=sum(self.bay_width)
        # Assigning values to instance variables
        self.no_of_bays = len(self.bay_width)-1            ######### number of bays
        self.no_of_stories = len(self.storey_height)-1     ######### number of stories
        self.no_of_elements_column = no_of_elements_column
        self.no_of_elements_beam = no_of_elements_beam
        self.no_of_nodes_column=self.no_of_elements_column-1
        self.no_of_nodes_beam=self.no_of_elements_beam-1
        self.node_to_node_height_column=[storey_height/self.no_of_elements_column for storey_height in self.storey_height]
        self.node_to_node_height_column=self.node_to_node_height_column[1:]
        self.node_to_node_length_beam=[bay_width/self.no_of_elements_beam for bay_width in self.bay_width]
        self.node_to_node_length_beam=self.node_to_node_length_beam[1:]
        self.beam_section=copy.deepcopy(beam_section)
        self.column_section=copy.deepcopy(column_section)   
        self.Main_Nodes=[]
        self.Leaning_Nodes=[]
        self.D_multiplier=load_combination_multipliers[0]      ### Dead Load multiplier
        self.L_multiplier=load_combination_multipliers[1]      ### Live Load multiplier
        self.L_r_multiplier=load_combination_multipliers[2]    ### Roof Live Load multiplier
        self.W_multiplier=load_combination_multipliers[3]      ### Wind Load multiplier
        self.Frame_id=Frame_id
        self.make_beam_section_detail_uniform()
        self.make_column_section_detail_uniform()
        self.load_timeseries_counter = 1
        self.load_pattern_counter = 1
        self.kwargs=kwargs

        defaults={'support':'All_Fixed',
                  'nip':4,
                  'mat_type':'Steel01',
                  'nfy':20,
                  'nfx':20,
                  'num_regions':10,
                  'D_floor_intensity':0,
                  'D_roof_intensity':0,
                  'L_floor_intensity':0,
                  'L_roof_intensity':0,
                  'Wind_load_floor':0,
                  'Wind_load_roof':0,
                  'Wall_load':0,
                  'Elastic_analysis':False,
                  'Second_order_effects':False,
                  'Notional_load':False,
                  'Geometric_Imperfection':False,
                  'Residual_Stress':True,
                  'stiffness_reduction':1,
                  'strength_reduction':1,
                  'geometric_imperfection_ratio':1/500,
                  'wind_load_dirn':'right',
                  'Leaning_column':True,
                  'Leaning_column_offset':4,
                  'Leaning_column_floor_load':0,
                  'Leaning_column_roof_load':0,
                  'plot_sections':False}
        
        for key,value in defaults.items():
            setattr(self,key,kwargs.get(key,value))
        
        if self.support=='All_Fixed':
            self.support_condition=['F']*(self.no_of_bays+1)
        else:
            self.support_condition=list(self.support)

        # Input validation
        if not (0<= self.no_of_bays < 10):
            raise ValueError("Number of bays must be from 1 to 9.")

        if not(1<=self.no_of_stories < 100):
            raise ValueError("Number of stories must be from 1 to 99.")

        if not(1<=self.no_of_elements_column < 10):
            raise ValueError("Number of elements in columns must be from 1 to 9.")

        if not(1<=self.no_of_elements_beam < 9):
            raise ValueError("Number of elements in beams must be from 1 to 9.")
        
        if len(self.support_condition)!= self.no_of_bays+1:
            raise ValueError(f"The number of arguments given for supports {self.support} should be equal to the number of columns, {(self.no_of_bays)+1}. ")
        

    def axis_i_floor_load(self, i):
        if i == 1:
            return [self.bay_i_floor_load(i)[0] * 0.5]
        elif i == (self.no_of_bays + 1):
            return [self.bay_i_floor_load(i - 1)[0] * 0.5]
        else:
            return [(self.bay_i_floor_load(i - 1)[0] + self.bay_i_floor_load(i)[0]) * 0.5]

    def axis_i_roof_load(self, i):
        if i == 1:
            return [self.bay_i_roof_load(i)[0] * 0.5]
        elif i == (self.no_of_bays + 1):
            return [self.bay_i_roof_load(i - 1)[0] * 0.5]
        else:
            return [(self.bay_i_roof_load(i - 1)[0] + self.bay_i_roof_load(i)[0]) * 0.5]

    def floor_notional_load(self):
        return (((self.D_multiplier*self.D_floor_intensity)+(self.L_multiplier*self.L_floor_intensity))*sum(self.bay_width))*0.002

    def roof_notional_load(self):
        return (((self.D_multiplier*self.D_roof_intensity)+(self.L_r_multiplier*self.L_roof_intensity))*sum(self.bay_width))*0.002


                # Control node for lateral deflection
    def get_control_node_and_dof(self,control_dir='L'):
        if control_dir=='L':
            control_direction='lateral'
            control_node = self.axis_i_roof_nodes(1)[0]   
            control_dof = 1  # horizontal displacement
        # Control node for vertical deflection of beam node
        else:
            control_direction='vertical'
            control_node = self.bay_i_internal_roof_nodes(1)[0]  
            control_dof = 2  # vertical displacement
        return control_node,control_dof 