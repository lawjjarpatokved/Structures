
import math as math
import pandas as pd
import openseespy.opensees as ops
import os
from libdenavit.OpenSees.plotting import plot_deformed_2d,plot_undeformed_2d,get_node_coords_and_disp
from libdenavit.OpenSees.get_fiber_data import *
from Plots import line_plot
from helpers import save_deformed_shape_frame,WF_Database
from math import pi, ceil
from libdenavit.section.wide_flange import *
from libdenavit.OpenSees import AnalysisResults
from libdenavit import find_limit_point_in_list, interpolate_list
import copy
import inspect
import opsvis
import matplotlib
matplotlib.use("TkAgg")
import numpy as np
import matplotlib.pyplot as plt
import opsvis as opsv
from pathlib import Path

#################################################
# density_of_steel=7850*kg/(m**3)
# g=9.81*m/(sec**2)
# E=29000*ksi
# G=77221*Mpa
# Fy=36*ksi
# Hk = 0.001*E           # Kinematic hardening modulus

# class Steel_Material:
#     def __init__(self,mat_tag,E,Fy,G,Hk,density):
#         self.mat_tag=mat_tag
#         self.E=E
#         self.Fy=Fy
#         self.G=G
#         self.Hk=Hk
#         self.density=density
#         self.b=Hk / (E + Hk) 

class Structures_2D:


    print_ops_status=True

    def __init__(self,width_of_bay,storey_height,
                no_of_elements_column, no_of_elements_beam,
                 beam_section,column_section,load_combination_multipliers,Frame_id,Material_obj,Steel_Grade,
                 **kwargs):
        # ---- Save a "constructor snapshot" for later cloning ---- This is useful for resetting or duplicating the model. Eg. Calculation of del2_over_del1
        self._init_spec = copy.deepcopy({k: v for k, v in locals().items() if k != "self"})
        # print(self._init_spec)
        # print("Constructor snapshot saved. You can use self._init_spec to clone or reset the model later.")
        # input()

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
        self.Material_obj=Material_obj
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
                  'Base_Wind_load':0,
                  'Wall_load':0,
                  'Elastic_analysis':False,
                  'Second_order_effects':False,
                  'Notional_load':False,
                  'Geometric_Imperfection':False,
                  'Residual_Stress':True,
                  'stiffness_reduction':1,
                  'strength_reduction':1,
                  'geometric_imperfection_ratio':1/500,
                  'initial_out_of_straightness_ratio':1/1000,
                  'wind_load_dirn':None,
                  'initial_out_of_straightness_dirn': None,
                  'Leaning_column':True,
                  'Leaning_column_offset':4,
                  'Leaning_column_floor_load':0,
                  'Leaning_column_roof_load':0,
                  'plot_sections':False,
                  'column_only_model':False,
                  'floor_nodes_free':False,
                  'wind_load_same_for_all_h':False,
                  'Steel_Grade':'Not_specified______??????????????'}
        
        for key,value in defaults.items():
            setattr(self,key,kwargs.get(key,value))


        self.make_beam_section_detail_uniform()

        self.make_column_section_detail_uniform()  

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
  
  

    @staticmethod
    def nth_digit(num, n, m=None):
        m=1 if m is None else m
        """
        Returns m digits starting from the nth digit (1-based index) of the number.
        For example, digits_from_n(123456, 2, 3) returns 234.
        """
        num_str = str(abs(num))  # Convert to string and ignore sign
        if n <= 0 or n > len(num_str):
            raise ValueError(f"n = {n} is out of bounds for number {num}")
        if n + m - 1 > len(num_str):
            raise ValueError(f"Requested {m} digits from position {n}, but number only has {len(num_str)} digits.")
        
        return int(num_str[n - 1:n - 1 + m])

    @staticmethod
    def find_bay_and_storey_for_beams(beam_connectivity):
        """
        Given beam connectivity as [tag, first_node, second_node],
        returns (bay, storey) based on encoded node number format.
        Assumes node number format encodes bay and storey as digits,
        e.g., node 203 means bay=2, storey=3.
        """
        first_node = beam_connectivity[0]
        bay = Structures_2D.nth_digit(first_node, 1)
        storey = Structures_2D.nth_digit(first_node, 2,2)
        # print(first_node,storey)
        return bay, storey
        
    @staticmethod
    def find_axis_and_storey_for_columns(column_connectivity):
        """
        Given beam connectivity as [tag, first_node, second_node],
        returns (bay, storey) based on encoded node number format.
        Assumes node number format encodes bay and storey as digits,
        e.g., node 203 means bay=2, storey=3.
        """
        first_node = column_connectivity[0]
        axis = Structures_2D.nth_digit(first_node, 1)
        storey = Structures_2D.nth_digit(first_node, 2,2) + 1  
        ''' 1 is added because my assumption is that column in storey
          1 means the column that support the 1st storey slab
        '''
        return axis, storey
    
    @staticmethod
    def generate_clean_csv_file_from_messy_out_files(path):
        # Read the file with all columns as string to avoid parsing errors
        df_raw = pd.read_csv(path, sep=r'\s+', header=None, dtype=str, engine='python')

        # Drop any empty rows
        df_raw.dropna(how='all', inplace=True)

        # Convert to numeric with errors coerced to NaN
        df = df_raw.apply(pd.to_numeric, errors='coerce')

        # Drop rows with any NaN (non-numeric values)
        df.dropna(inplace=True)

        if df.empty:
            print(f" No valid numeric data found in: {path}")
            return

        # Select the last row and convert to inches
        df_clean = (df.iloc[-1] / 0.0254).to_frame()

        # Reset index and label the column
        df_clean.reset_index(drop=True, inplace=True)
        df_clean.columns = ['Value_in_inches']

        # Save as clean CSV
        folder = os.path.dirname(path)
        base_name = os.path.splitext(os.path.basename(path))[0]
        output_path = os.path.join(folder, base_name + '_clean.csv')
        df_clean.to_csv(output_path, index=False)

        print(f" Saved clean displacement to: {output_path}")


    def generate_Nodes_and_Element_Connectivity(self):
        x_coord=0
        for i in range (self.no_of_bays+1): ##### 0 to 8
            x_coord= x_coord+self.bay_width[i]
            y_coord=0
            for j in range(self.no_of_stories+1): # 0 to 98
                y_coord= y_coord+ self.storey_height[j]
                node_tag=(i+1)*100+j
                self.Main_Nodes.append([node_tag,x_coord,y_coord])

        self.NODES_TO_FIX= [Nodes for Nodes in self.Main_Nodes if Nodes[2]==0]
        self.main_nodes_except_base=[Nodes for Nodes in self.Main_Nodes if Nodes[2]!=0]

        

########### This code chunk generates the additional nodes and connectivity of columns based on the number of elements required ############
        temp_column_node_pairs = []

        for node1 in self.Main_Nodes:
            tag1 = node1[0]
            x1=node1[1]
            y1=node1[2]
            for node2 in  self.Main_Nodes:
                tag2 = node2[0]
                x2=node2[1]
                y2=node2[2]
                if abs(tag1 - tag2) == 1 and (tag1 // 100 == tag2 // 100) and tag2>tag1:
                    pair = [tag1,x1,y1, tag2]  # ensures lower, upper
                    if pair not in temp_column_node_pairs:
                        temp_column_node_pairs.append(pair)

        self.column_intermediate_nodes=[]
        self.column_connectivity=[]
        element_tag=1
        for node in temp_column_node_pairs:
            temp_column_connectivity=[node[0]]
            for i in range(self.no_of_nodes_column):
                x = node[1]
                y = node[2] + (i+1) * self.node_to_node_height_column[node[0]%100]
                tag = node[0] * 10 + i
                self.column_intermediate_nodes.append([tag, x, y])
                temp_column_connectivity.append(tag)
            temp_column_connectivity.append(node[3])
            for i in range(len(temp_column_connectivity)-1):
                self.column_connectivity.append([element_tag,temp_column_connectivity[i],temp_column_connectivity[i+1]])
                element_tag+=1

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################



        temp_beam_node_pairs = []

        for node1 in self.Main_Nodes:
            tag1 = node1[0]
            x1=node1[1]
            y1=node1[2]
            for node2 in  self.Main_Nodes:
                tag2 = node2[0]
                x2=node2[1]
                y2=node2[2]
                if abs(tag1 - tag2) == 100 and y1!=0 and y2!=0 and tag2>tag1:
                    pair = [tag1,x1,y1, tag2]  # ensures lower, upper
                    if pair not in temp_beam_node_pairs:
                        temp_beam_node_pairs.append(pair)


        self.beam_intermediate_nodes=[]
        self.beam_connectivity=[]

        for node in temp_beam_node_pairs:
            temp_beam_connectivity=[node[0]]
            for i in range(self.no_of_nodes_beam):
                x = node[1]+ (i+1) * self.node_to_node_length_beam[(node[0]//100)-1]
                y = node[2] 
                tag = node[0] * 100 + i
                self.beam_intermediate_nodes.append([tag, x, y])
                temp_beam_connectivity.append(tag)
            temp_beam_connectivity.append(node[3])
            for i in range(len(temp_beam_connectivity)-1):
                self.beam_connectivity.append([element_tag,temp_beam_connectivity[i],temp_beam_connectivity[i+1]])
                element_tag+=1
       
        self.beam_nodes = [node for node in self.Main_Nodes if node not in self.NODES_TO_FIX]
        self.beam_nodes += self.beam_intermediate_nodes
        
        self.all_nodes=self.Main_Nodes+self.column_intermediate_nodes+self.beam_intermediate_nodes


        #### Generating leaning column nodes and connectivity if required ####
        if self.Leaning_column:
            x_coord=self.length_of_frame + self.Leaning_column_offset
            y_coord=0
            for j in range(self.no_of_stories+1): # 0 to 98
                y_coord= y_coord+ self.storey_height[j]
                node_tag=100000+j
                self.Leaning_Nodes.append([node_tag,x_coord,y_coord])
            

            ## Generate nodal connectivity for leaning columns
            self.leaning_column_connectivity=[]
            ## assigning element tags continuing from last element tag
            for i in range(self.no_of_stories):
                element_tag+=1
                self.leaning_column_connectivity.append([element_tag,self.Leaning_Nodes[i][0],self.Leaning_Nodes[i+1][0]])

            ## Generate nodal connectivity between main frame nodes and leaning column nodes at each floor level to later use in defining constraints
            self.leaning_column_and_frame_beam_connectivity=[]
            for i in range(1,self.no_of_stories+1):
                element_tag+=1
                main_frame_node_tag=(self.no_of_bays+1)*100 + i
                self.leaning_column_and_frame_beam_connectivity.append([element_tag,main_frame_node_tag,self.Leaning_Nodes[i][0]])


###################### The code below is written so as to change the order of elements for direct comparison with Ziemian results
        self.sorted_column_connectivity= sorted(
                        self.column_connectivity,key=lambda col: tuple(reversed(Structures_2D.find_axis_and_storey_for_columns(col[1:])))
                        )
        
        self.column_member_list = []
        group_size=self.no_of_elements_column
        for i in range(0, len(self.sorted_column_connectivity), group_size):
            group = self.sorted_column_connectivity[i:i+group_size]
            if len(group) == group_size:
                member_tag = (i // group_size) + 1
                member_element_tags = [e[0] for e in group]
                self.column_member_list.append([member_tag] + member_element_tags)


        self.sorted_beam_connectivity=sorted(self.beam_connectivity,key=lambda b: Structures_2D.find_bay_and_storey_for_beams(b[1:])[1])

        self.beam_member_list = []
        group_size = self.no_of_elements_beam  

        for i in range(0, len(self.sorted_beam_connectivity), group_size):
            group = self.sorted_beam_connectivity[i:i+group_size]
            if len(group) == group_size:
                member_tag = len(self.column_member_list) + (i // group_size) + 1  # Continues from last column member
                member_element_tags = [e[0] for e in group]
                self.beam_member_list.append([member_tag] + member_element_tags)

        self.sorted_element_connectivity=self.sorted_column_connectivity+self.sorted_beam_connectivity
        self.member_list=self.column_member_list+self.beam_member_list

    def bay_i_internal_floor_nodes(self, i):
        '''This function returns the node tags of all the floor nodes that lie in ith bay.This is only applicable 
          to Moment Frames not column models  '''
        return [
            node[0]
            for node in self.beam_nodes
            if len(str(node[0])) == 5
            and Structures_2D.nth_digit(node[0], 1) == i
            and Structures_2D.nth_digit(node[0], 2,2) != self.no_of_stories
        ]
    
    def bay_i_internal_roof_nodes(self, i):
        '''This function returns the node tags of all the roof nodes that lie in ith bay.This is only applicable 
          to Moment Frames not column models  '''
        return [
            node[0]
            for node in self.beam_nodes
            if len(str(node[0])) == 5
            and Structures_2D.nth_digit(node[0], 1) == i
            and Structures_2D.nth_digit(node[0], 2,2) == self.no_of_stories
        ]

    def axis_i_roof_nodes(self, i):
        '''This function returns the node tags of all the roof nodes that lie in ith axis.Applicable to both Moment Frames and Column models.
         For columns, set i=1 as there is only one axis'''
        return [
            node[0]
            for node in self.beam_nodes
            if len(str(node[0])) == 3
            and Structures_2D.nth_digit(node[0], 1) == i
            and Structures_2D.nth_digit(node[0], 2,2) == self.no_of_stories
        ]
    
    def axis_i_floor_nodes(self, i):
        '''This function returns the node tags of all the floor nodes that lie in ith axis.Applicable to both Moment Frames and Column models.
         For columns, set i=1 as there is only one axis  '''
        return [
            node[0]
            for node in self.beam_nodes
            if len(str(node[0])) == 3
            and Structures_2D.nth_digit(node[0], 1) == i
            and Structures_2D.nth_digit(node[0], 2,2) != self.no_of_stories
        ]

    def bay_i_floor_load(self, i):
        ''' This function returns the point loads to be applied to all the internal floor nodes of ith bay.This is only applicable 
          to Moment Frames not column models  '''
        return [((-self.D_multiplier*self.D_floor_intensity)+(-self.L_multiplier*self.L_floor_intensity))*(self.bay_width[i]/self.no_of_elements_beam)]
    
    def bay_i_roof_load(self, i):
        ''' This function returns the point loads to be applied to all the internal roof nodes of ith bay.This is only applicable 
          to Moment Frames not column models  '''
        return [((-self.D_multiplier*self.D_roof_intensity)+(-self.L_r_multiplier*self.L_roof_intensity))*(self.bay_width[i]/self.no_of_elements_beam)]

    def create_distorted_nodes_and_element_connectivity(self,geometric_imperfection_ratio=None,initial_out_of_straightness_ratio=None):
        if self.Geometric_Imperfection:
            print('Working in imperfect geometry')
            ratio = (geometric_imperfection_ratio if geometric_imperfection_ratio is not None else self.geometric_imperfection_ratio) * (1 if self.wind_load_dirn=='right' else -1)
            for i in range(len(self.all_nodes)):
                    self.all_nodes[i][1] = self.all_nodes[i][1] +ratio * self.all_nodes[i][2]


            ## Initial out of straightness
            for members in self.column_member_list:
                member_tag=members[0]
                member_elements=members[1:]
                member_connectivity = [connectivity for connectivity in self.sorted_column_connectivity if connectivity[0] in member_elements]
                member_nodes = [member_connectivity[0][1]] + [connectivity[2] for connectivity in member_connectivity]
                node_dict = {node[0]: node for node in self.all_nodes}
                # Get member end coordinates
                x1, y1 = node_dict[member_nodes[0]][1], node_dict[member_nodes[0]][2]
                x2, y2 = node_dict[member_nodes[-1]][1], node_dict[member_nodes[-1]][2]

                # Member length
                L = np.hypot(x2 - x1, y2 - y1) 

                for node_tag in member_nodes:
                    node = node_dict[node_tag]

                    # Distance of node along the member
                    s = np.hypot(node[1] - x1, node[2] - y1)

                    # Initial out-of-straightness
                    initial_out_of_straightness=(initial_out_of_straightness_ratio if initial_out_of_straightness_ratio is not None else self.initial_out_of_straightness_ratio) * (1 if self.initial_out_of_straightness_dirn=='right' else -1)
                    imperfection = initial_out_of_straightness*L * np.sin(np.pi * s / L)

                    # Modify x-coordinate directly in self.all_nodes
                    node[1] += imperfection
    
        else:
            print('Working in nominal geometry')
            pass

    def make_beam_section_detail_uniform(self):
        self.beam_case = list(self.beam_section.keys())[0]
        nested_dict = self.beam_section[self.beam_case]

        self.no_of_beam_sections = len(nested_dict)

        # Assign unique tags to unique section names
        unique_sections = sorted(set(nested_dict.values()))
        self.beam_section_tags = {section: idx + 1 for idx, section in enumerate(unique_sections)}

        # Update the nested dictionary to include the tag
        updated_nested_dict = {
            key: (value, self.beam_section_tags[value]) for key, value in nested_dict.items()
        }

        # Update the main beam_section
        self.beam_section[self.beam_case] = updated_nested_dict

    def make_column_section_detail_uniform(self):
        self.column_case = list(self.column_section.keys())[0]
        nested_dict = self.column_section[self.column_case]

        self.no_of_column_sections = len(nested_dict)

        # Extract unique section names (still keep axis handling separately)
        unique_sections = sorted(set(value[0] for value in nested_dict.values()))

        # Assign tag numbers (tags only for now)
        section_tags = {
            section: idx + 1 + self.no_of_beam_sections
            for idx, section in enumerate(unique_sections)
        }

        # Initialize updated dicts
        updated_nested_dict = {}
        self.column_section_tags = {}

        # Go through each column and assign (section, tag, axis)
        
        for key, (section, axis) in nested_dict.items():
            tag = section_tags[section]
            updated_nested_dict[key] = (section, tag, axis)

            # Store axis info alongside tag — update only if not already set
            # or overwrite to latest axis if duplicates exist
            self.column_section_tags[section] = (tag, axis)

        self.column_section[self.column_case] = updated_nested_dict

    def build_ops_model(self):
        ops.wipe()
        
        ops.model('basic','-ndm',2,'-ndf',3)
        
        for single_node in self.all_nodes:
            ops.node(single_node[0],single_node[1],single_node[2])
    
        # ops.printModel()

        for base_nodes, support_condition in zip(self.NODES_TO_FIX, self.support_condition):
            node_tag = base_nodes[0]
            if support_condition.upper() == 'F':       # Fixed: UX, UY, RZ all fixed
                ops.fix(node_tag, 1, 1, 1)
            elif support_condition.upper() == 'P':     # Pinned: UX and UY fixed, RZ free
                ops.fix(node_tag, 1, 1, 0)
            else:                                      # Wrong condition
                raise ValueError(f"Unsupported support condition{support_condition}.Expected 'F' or 'P'.")
            
        if self.column_only_model and not self.floor_nodes_free:
            for nodes in self.main_nodes_except_base:
                node_tag=nodes[0]
                ops.fix(node_tag,0,0,1)    

        if self.Leaning_column:
            ## Define leaning column nodes
            for leaning_column_node in self.Leaning_Nodes:
                ops.node(leaning_column_node[0],leaning_column_node[1],leaning_column_node[2])

            ## Fix the rotational dof of leaning column nodes except the base node
            for leaning_node in self.Leaning_Nodes[1:]:
                ops.fix(leaning_node[0],0,0,1)          # Restricting only rotational dof for leaning column nodes except base node

            ## Pin the base of leaning column if it exists  
            base_leaning_node_tag=self.Leaning_Nodes[0][0]
            ops.fix(base_leaning_node_tag,1,1,1)          # Pinned support for leaning column base node


        col_and_beam_TransTag = 1

        if self.Elastic_analysis:
            mat_type = 'Elastic'
            frc = 0
        else:
            mat_type = self.mat_type
            frc = -0.3 * self.Material_obj.Fy if self.Residual_Stress else 0

    

        # Beams
        for beam_section_name, beam_section_tag in self.beam_section_tags.items():
            beam_data = WF_Database(beam_section_name)
            beam = I_shape(beam_data.d, beam_data.tw, beam_data.bf, beam_data.tf,   
                        Fy=self.Material_obj.Fy, E=self.Material_obj.E,
                        A=beam_data.A, 
                        Ix=beam_data.Ix,Zx=beam_data.Zx,Sx=beam_data.Sx,rx=beam_data.rx,
                        Iy=beam_data.Iy,Zy=beam_data.Zy,Sy=beam_data.Sy,ry=beam_data.ry,
                        J=beam_data.J,Cw=beam_data.Cw,rts=beam_data.rts,ho=beam_data.ho)
            beam.build_ops_fiber_section(beam_section_tag,
                                        start_material_id=self.Material_obj.mat_tag,
                                        mat_type=mat_type,
                                        nfy=self.nfy, nfx=self.nfx,
                                        frc=frc,num_regions=self.num_regions,
                                        stiffness_reduction=self.stiffness_reduction,strength_reduction=self.strength_reduction,
                                        axis='x')
            self.Material_obj.mat_tag += 2 * self.num_regions + 2
            ops.beamIntegration("Lobatto", beam_section_tag, beam_section_tag, self.nip)
            setattr(self, beam_section_name, beam)

        # Columns
        for column_section_name, (column_section_tag, axis) in self.column_section_tags.items():
            column_data = WF_Database(column_section_name)
            column = I_shape(column_data.d, column_data.tw, column_data.bf, column_data.tf,
                            Fy=self.Material_obj.Fy, E=self.Material_obj.E,
                            A=column_data.A,
                            Ix=column_data.Ix, Zx=column_data.Zx, Sx=column_data.Sx, rx=column_data.rx,
                            Iy=column_data.Iy, Zy=column_data.Zy, Sy=column_data.Sy, ry=column_data.ry,
                            J=column_data.J, Cw=column_data.Cw, rts=column_data.rts, ho=column_data.ho)
            column.build_ops_fiber_section(column_section_tag,
                                        start_material_id=self.Material_obj.mat_tag,
                                        mat_type=mat_type,
                                        nfy=self.nfy, nfx=self.nfx,
                                        frc=frc,num_regions=self.num_regions,
                                        stiffness_reduction=self.stiffness_reduction,strength_reduction=self.strength_reduction,
                                        axis=axis)
            self.Material_obj.mat_tag += 2 * self.num_regions + 2
            ops.beamIntegration("Lobatto", column_section_tag, column_section_tag, self.nip)
            setattr(self, column_section_name, column)


        ########### check whether Second order effects need to be included or not ############
        if self.Second_order_effects:    
            ops.geomTransf("PDelta", col_and_beam_TransTag)

        else:
            ops.geomTransf("Linear", col_and_beam_TransTag)
 





        # ----- Columns -----
        for i, column_ij_node in enumerate(self.column_connectivity):
            
            axis, storey = Structures_2D.find_axis_and_storey_for_columns(column_ij_node[1:])

            if self.column_case == 'common_and_exceptions':
                key = f'({axis},{storey})'
                if key in self.column_section[self.column_case]:
                    # print(f"node {column_ij_node[1]} and {column_ij_node[2]} are in axis {axis} and storey {storey}, which is an exception.")
                    section_name = self.column_section[self.column_case][key][0]
                    section_tag  = self.column_section[self.column_case][key][1]
                    bending_axis = self.column_section[self.column_case][key][2]
                else:
                    # print(f"node {column_ij_node[1]} and {column_ij_node[2]} are not in axis {axis} and storey {storey}, which is common.")
                    section_name = self.column_section[self.column_case]['common'][0]
                    section_tag  = self.column_section[self.column_case]['common'][1]
                    bending_axis = self.column_section[self.column_case]['common'][2]

            elif self.column_case == 'same_for_storey':
                key = str(storey)
                if key in self.column_section[self.column_case]:
                    section_name = self.column_section[self.column_case][key][0]
                    section_tag  = self.column_section[self.column_case][key][1]
                    bending_axis = self.column_section[self.column_case][key][2]
                else:
                    raise KeyError(f"No column section defined for storey {storey}")

            else:
                raise ValueError(f"Unsupported column_case: {self.column_case}. Possible issues: lowercase, spelling error.")

            # Create column element
            eleTag = column_ij_node[0]
            ops.element('forceBeamColumn', eleTag, *column_ij_node[1:3], col_and_beam_TransTag, section_tag,'-mass', 1)

            # Append section name to column entry
            self.column_connectivity[i] = column_ij_node + [section_name] +[bending_axis]+['col']





        # ----- Beams -----
        for i, beam_ij_node in enumerate(self.beam_connectivity):
            bay, storey = Structures_2D.find_bay_and_storey_for_beams(beam_ij_node[1:])

            if self.beam_case == 'common_and_exceptions':
                key = f'({bay},{storey})'
                if key in self.beam_section[self.beam_case]:
                    # print(f"node {beam_ij_node[1]} and {beam_ij_node[2]} are in bay {bay} and storey {storey}, which is an exception.")
                    section_name = self.beam_section[self.beam_case][key][0]
                    section_tag  = self.beam_section[self.beam_case][key][1]
                else:
                    # print(f"node {beam_ij_node[1]} and {beam_ij_node[2]} are not in bay {bay} and storey {storey}, which is common.")
                    section_name = self.beam_section[self.beam_case]['common'][0]
                    section_tag  = self.beam_section[self.beam_case]['common'][1]

            elif self.beam_case == 'same_for_storey':
                key = str(storey)
                if key in self.beam_section[self.beam_case]:
                    section_name = self.beam_section[self.beam_case][key][0]
                    section_tag  = self.beam_section[self.beam_case][key][1]
                else:
                    raise KeyError(f"No beam section defined for storey {storey}")

            else:
                raise ValueError(f"Unsupported beam_case: {self.beam_case}, Possible Errors: Lower case, spelling")

            # Create beam element
            eleTag = beam_ij_node[0]
            ops.element('forceBeamColumn', eleTag, *beam_ij_node[1:3], col_and_beam_TransTag, section_tag,'-mass', 1)
            # print(f"Defined element {eleTag} between nodes {column_ij_node[1]} and {column_ij_node[2]}")

            # Append section name to beam entry
            self.beam_connectivity[i] = beam_ij_node + [section_name] + ['x']+['beam']   ###'x' is hard coded because the beams are expected to be bending about major axis only (Discussed in meeting with prof.)

        self.all_element_connectivity_section_and_bending_axes_detail=self.column_connectivity+self.beam_connectivity

        ## Define leaning column elements (corrotational truss/truss with large AE/L) if required ##
        if self.Leaning_column:
            for leaning_col_ij in self.leaning_column_connectivity:
                eleTag = leaning_col_ij[0]
                node_i = leaning_col_ij[1]
                node_j = leaning_col_ij[2]
                if self.Second_order_effects:
                    ops.element('corotTruss', eleTag, node_i, node_j, 100 ,1,'-rho',1)
                else:
                    ops.element('truss', eleTag, node_i, node_j, 100 ,1, '-rho',1)

            ## Define geometric constraint between leaning column nodes and main frame nodes at each floor level ## 
            for constraint in self.leaning_column_and_frame_beam_connectivity:
                eleTag = constraint[0]
                main_frame_node_tag = constraint[1]
                leaning_column_node_tag = constraint[2]
                ops.equalDOF(main_frame_node_tag,leaning_column_node_tag,1)  # Constraining UX and UY of leaning column node to main frame node




        self.roof_beams=[]
        for beams in self.beam_connectivity:
            if Structures_2D.nth_digit(beams[1],2,2)==self.no_of_stories:
                self.roof_beams.append(beams)

    def rebuild_with_overrides(self, message=False, **overrides ):

        if message:
            print("WARNING: You are using the rebuild_with_overrides method, which creates a new instance of the class with overridden parameters." \
            " Make sure to assign the result to a variable, e.g., new_instance = old_instance.rebuild_with_overrides(**overrides). Furthermore, if you are using this" \
            "method for the first time, make sure that you have given correct override parameters. You can do a simple chekck by comparing the output of the original instance and the new instance for a simple case where you know the expected output. ")

            input("Press Enter to continue...")

        ## This function rebuilds the current instance with overridden parameters.This is 
        ## useful for parametric studies where only a few parameters need to be changed.
        ## For example, this method can be used for the calculation of del2_over_del1 parameter
        ## for a frame whatever be the original analysis choice (Second order effects: True or False) made during instantiation of the class.
        spec = copy.deepcopy(self._init_spec)
        spec.setdefault("kwargs", {})

        sig = inspect.signature(self.__class__.__init__)
        explicit_params = {
            p.name for p in sig.parameters.values()
            if p.name not in ("self",) and p.kind != inspect.Parameter.VAR_KEYWORD
        }
        # print("Explicit parameters:", explicit_params)
    
        # Apply overrides: explicit args stay explicit, everything else goes into **kwargs

        # print('Overrides to apply:', overrides)
        for k, v in overrides.items():
            if k in explicit_params:
                spec[k] = copy.deepcopy(v)
            else:
                spec["kwargs"][k] = copy.deepcopy(v)

        #  split the call into explicit args + **kwargs dict
        kwargs_dict = spec.pop("kwargs", {})
        # print("Spec after overrides:", spec)
        # print("Kwargs dict:", kwargs_dict)
        
        return self.__class__(**spec, **kwargs_dict)
    
    def return_max_of_fiber_strain_in_all_elements(self):
        ## returns the maximum strain from among all elements in the Frame
        maximum_compression_strain=[]
        maximum_tensile_strain=[]
        for ele in self.all_element_connectivity_section_and_bending_axes_detail:
            ele_tag=ele[0]
            # ele_node_i=ele[1]
            # ele_node_j=ele[2]
            ele_section_name=ele[3]        #### str  'W8X15'
            ele_bending_axis=ele[4]        #### str  'x'
            section_obj = getattr(self, ele_section_name)   # retrieves the I_shape instance
            d = section_obj.d
            bf = section_obj.bf

            compression_strain = []
            tensile_strain = []

            for i in range(self.nip):
                axial_strain, curvatureX, curvatureY = 0, 0, 0
                if ele_bending_axis=='x':
                    axial_strain, curvatureX = ops.eleResponse(ele_tag,  # element tag
                                                                'section', i+1,  # select integration point
                                                                'deformation')  # response type               
                elif ele_bending_axis=='y':
                    axial_strain, curvatureY = ops.eleResponse(ele_tag,  # element tag
                                                                'section', i+1,  # select integration point
                                                                'deformation')  # response type
                else:
                    raise ValueError("The axis is not supported.")
                
                compression_strain.append(section_obj.maximum_compression_strain(axial_strain,curvatureX,curvatureY))
                tensile_strain.append(section_obj.maximum_tensile_strain(axial_strain,curvatureX,curvatureY))
            # print('Element tag',ele_tag)
            # print('Section_name',ele_section_name, ele_bending_axis)
            # print('d',d)
            # print('bf',bf)
            # print('Axial_strain',axial_strain)
            # print('CurvatureX',curvatureX)
            # print('CurvatureY',curvatureY)
            # print(compression_strain)   
            # print(tensile_strain)
            maximum_compression_strain.append(max(compression_strain,key=lambda x:abs(x)))
            # print(maximum_compression_strain)
            maximum_tensile_strain.append(max(tensile_strain,key=lambda x:abs(x)))
            # print(maximum_tensile_strain)
        # print('Maximum among all elements')
        # print("642",maximum_compression_strain)
        # print('643',maximum_tensile_strain)
        max_abs_compression_strain= max(abs(c) for c in maximum_compression_strain)
        max_abs_tensile_strain= max(abs(t) for t in maximum_tensile_strain)
        # print(max_abs_compression_strain)
        # print(max_abs_tensile_strain)

        return max(max_abs_compression_strain,max_abs_tensile_strain)

    def return_P_M_M_interaction_values(self):
        P_M_M_interaction_all_elements = []
        Element_Forces=[]
        for ele in self.all_element_connectivity_section_and_bending_axes_detail:
            ele_tag = ele[0]
            ele_node_i = ele[1]
            ele_node_j = ele[2]
            ele_section_name = ele[3]        # e.g., 'W8X15'
            ele_bending_axis = ele[4]        # e.g., 'x'
            ele_type = ele[5]

            coords_i = ops.nodeCoord(ele_node_i)
            coords_j = ops.nodeCoord(ele_node_j)
            L = math.sqrt((coords_j[0] - coords_i[0])**2 + (coords_j[1] - coords_i[1])**2)

            section_obj = getattr(self, ele_section_name)
            member_name = f"member_{ele_tag}_{ele_section_name}_{ele_bending_axis}"

            member_obj = WideFlangeMember_AISC2022(
                section_obj,
                Fy=section_obj.Fy,
                E=section_obj.E,
                L=L
            )
            setattr(self, member_name, member_obj)

            # Element forces
            forces = ops.eleResponse(ele_tag, 'localForce')
            Pr = abs(forces[0])
            Mr_i = forces[2]
            Mr_j = forces[5]
            max_Mr = max(abs(Mr_i), abs(Mr_j))

            Mrx = Mry = 0
            Mcx = member_obj.Mnx(Lb=0, Cb=1)
            Mcy = member_obj.Mny()

            if ele_bending_axis == 'x':
                Mrx = max_Mr
                Leff = L * (self.no_of_elements_column if ele_type == 'col' else self.no_of_elements_beam)
                Pcc = member_obj.Pnc(Lcx=Leff, Lcy=0)

            elif ele_bending_axis == 'y':
                Mry = max_Mr
                Leff = L * (self.no_of_elements_column if ele_type == 'col' else self.no_of_elements_beam)
                Pcc = member_obj.Pnc(Lcx=0, Lcy=Leff)

            # Interaction Eqn H1-1a or H1-1b
            if Pr / Pcc >= 0.2:
                P_M_M_interaction = Pr / Pcc + (8 / 9) * ((Mrx / Mcx) + (Mry / Mcy))
            else:
                P_M_M_interaction = Pr / (2 * Pcc) + ((Mrx / Mcx) + (Mry / Mcy))

            # Store as tuple
            P_M_M_interaction_all_elements.append((ele_tag, P_M_M_interaction))
            Element_Forces.append((ele_tag,forces))

        # Find max interaction and its element tag
        max_ele_tag, max_PMM = max(P_M_M_interaction_all_elements, key=lambda x: x[1])

        return max_PMM, max_ele_tag,P_M_M_interaction_all_elements,Element_Forces

    def add_vertical_dead_live_wall_notional_loads(self,vertical_load_scale=1):   
        """
        Adds dead, live, wall loads to the model, scaled by a user-defined load_scale factor.
        
        Parameters:
            load_scale (float): Scaling factor for all loads. Use 1.0 for full load, <1.0 for partial.
            This is helpful when performing load controlled or displacement controlled analysis. While 
            performing displacement controlled analysis, it is better to apply small portion of the load
            so that the load factors are positive throughout the analysis.
        """

        ##Vertical Loads
        # Node-based Dead and Live Loads 
        for bay in range(1, self.no_of_bays + 1):
            loaded_nodes_floor = self.bay_i_internal_floor_nodes(i=bay)
            load_value_floor = self.bay_i_floor_load(i=bay)
            for node in loaded_nodes_floor:
                ops.load(node, 0.0, vertical_load_scale * load_value_floor[0], 0.0)

            loaded_nodes_roof = self.bay_i_internal_roof_nodes(bay)
            load_value_roof = self.bay_i_roof_load(bay)
            for node in loaded_nodes_roof:
                ops.load(node, 0.0, vertical_load_scale * load_value_roof[0], 0.0)

        for axis in range(1, self.no_of_bays + 2):
            loaded_nodes_floor = self.axis_i_floor_nodes(axis)
            load_value_floor = self.axis_i_floor_load(axis)
            for node in loaded_nodes_floor:
                ops.load(node, 0.0, vertical_load_scale * load_value_floor[0], 0.0)

            loaded_nodes_roof = self.axis_i_roof_nodes(axis)
            load_value_roof = self.axis_i_roof_load(axis)
            for node in loaded_nodes_roof:
                ops.load(node, 0.0, vertical_load_scale * load_value_roof[0], 0.0)

        # Wall Load
        for node in self.axis_i_floor_nodes(1):
            wall_load = self.Wall_load * self.D_multiplier
            ops.load(node, 0, -vertical_load_scale * wall_load, 0)

        for node in self.axis_i_floor_nodes(self.no_of_bays + 1):
            wall_load = self.Wall_load * self.D_multiplier
            ops.load(node, 0, -vertical_load_scale * wall_load, 0)

        for node in self.axis_i_roof_nodes(1):
            wall_load = self.Wall_load * self.D_multiplier
            ops.load(node, 0, -vertical_load_scale * wall_load / 2, 0)

        for node in self.axis_i_roof_nodes(self.no_of_bays + 1):
            wall_load = self.Wall_load * self.D_multiplier
            ops.load(node, 0, -vertical_load_scale * wall_load / 2, 0)


        if self.Notional_load:
            if self.wind_load_dirn.lower() == 'right':
                for node in self.axis_i_floor_nodes(1):
                    ops.load(node, vertical_load_scale * self.floor_notional_load(), 0, 0.0)
                    # print('line852',node,self.floor_notional_load())
                for node in self.axis_i_roof_nodes(1):
                    ops.load(node, vertical_load_scale * self.roof_notional_load(), 0, 0.0)
                    # print('line856',node,lateral_load_scale * self.roof_notional_load()) 

            elif self.wind_load_dirn.lower() == 'left':
                for node in self.axis_i_floor_nodes(self.no_of_bays + 1):
                    ops.load(node, -vertical_load_scale * self.floor_notional_load(), 0, 0.0)
                    # print('line862',node,self.floor_notional_load())
                for node in self.axis_i_roof_nodes(self.no_of_bays + 1):
                    ops.load(node, -vertical_load_scale * self.roof_notional_load(), 0, 0.0)
                    # print('line866',node,lateral_load_scale * self.roof_notional_load())

        if self.Leaning_column:
            ## Add vertical loads to the floor and roof nodes of leaning columns ##
            ## floor nodes ##
            for i in range(1,self.no_of_stories):
                ops.load(self.Leaning_Nodes[i][0], 0.0, -vertical_load_scale * self.Leaning_column_floor_load, 0.0)
            ## roof node ##
            ops.load(self.Leaning_Nodes[-1][0], 0.0, -vertical_load_scale * self.Leaning_column_roof_load, 0.0)

    def add_lateral_wind_loads(self, lateral_load_scale=1.0):
        ##Lateral Loads
        # Wind  Load 
        
        if self.wind_load_dirn is None:
            print('Lateral Loads not applied in the model')

        elif self.wind_load_dirn.lower() == 'right':
            print('Applying wind loads towards right')

            loaded_nodes=self.axis_i_floor_nodes(1)+self.axis_i_roof_nodes(1)
            loaded_nodes.sort()

            if self.wind_load_same_for_all_h:
                for node in loaded_nodes:
                    ops.load(node, lateral_load_scale * self.Base_Wind_load * self.W_multiplier, 0, 0.0)

            else:
                for node in loaded_nodes:
                    coords = ops.nodeCoord(node)
                    H = ops.nodeCoord(loaded_nodes[0])[1]
                    wind_load=self.Base_Wind_load*coords[1]/H
                    ops.load(node, lateral_load_scale * wind_load * self.W_multiplier, 0, 0.0)


        elif self.wind_load_dirn.lower() == 'left':
            print('Applying wind loads towards left')

            loaded_nodes=self.axis_i_floor_nodes(self.no_of_bays + 1)+self.axis_i_roof_nodes(self.no_of_bays + 1)
            loaded_nodes.sort()

            if self.wind_load_same_for_all_h:
                for node in loaded_nodes:
                    ops.load(node, -lateral_load_scale * self.Base_Wind_load* self.W_multiplier, 0, 0.0)
            else:
                for node in loaded_nodes:
                    coords = ops.nodeCoord(node)
                    H = ops.nodeCoord(loaded_nodes[0])[1]
                    wind_load=self.Base_Wind_load*coords[1]/H
                    ops.load(node, -lateral_load_scale *wind_load * self.W_multiplier, 0, 0.0)




        
                
        # Visualization disabled to avoid GUI issues in non-interactive environments
        # opsv.plot_model()
        # opsv.plot_load()

    def return_drift_of_all_storeys_at_given_axis(self,i): 
        ## Calculate the drift of each storey at one of the axis (say axis i).
        ## At each axis I need to separately find the drift of first story using the nodes to fix and 
        ## axis i floor nodes. and the use loop to find the drift of other storeys. and again find the drift
        ## of roof using axis i roof nodes and last element of axis i floor nodes.
        drift = []

        base_node_tag = self.NODES_TO_FIX[i-1][0]
        storey_nodes_tags = self.axis_i_floor_nodes(i)
        roof_node_tag = self.axis_i_roof_nodes(i)[0]

        # -------- Special Case: One Storey --------
        if self.no_of_stories == 1:
            single_storey_drift = (
                ops.nodeDisp(roof_node_tag, 1)
                - ops.nodeDisp(base_node_tag, 1)
            )
            drift.append(single_storey_drift)
            return drift
        # ------------------------------------------

        ## First storey drift
        first_storey_drift = (
            ops.nodeDisp(storey_nodes_tags[0], 1)
            - ops.nodeDisp(base_node_tag, 1)
        )
        drift.append(first_storey_drift)

        ## Other storeys drift
        for storey in range(1, self.no_of_stories - 1):
            lower_storey_node_tag = storey_nodes_tags[storey - 1]
            upper_storey_node_tag = storey_nodes_tags[storey]
            storey_drift = (
                ops.nodeDisp(upper_storey_node_tag, 1)
                - ops.nodeDisp(lower_storey_node_tag, 1)
            )
            drift.append(storey_drift)

        ## Roof drift
        last_storey_node_tag = storey_nodes_tags[-1]
        roof_drift = (
            ops.nodeDisp(roof_node_tag, 1)
            - ops.nodeDisp(last_storey_node_tag, 1)
        )
        drift.append(roof_drift)
        print(drift)
        return drift       

    def run_load_controlled_anlaysis(self,**kwargs):

        

        incr_LCA= kwargs.get('incr_LCA', 0.1)           ######### LCA refers to Load Controlled Analysis
        num_steps_LCA= kwargs.get('num_steps_LCA', 200)            ######### LCA refers to Load Controlled Analysis
        steel_strain_limit = kwargs.get('steel_strain_limit', 0.05)
        eigenvalue_limit = kwargs.get('eigenvalue_limit', 0)
        P_M_M_interaction_limit=kwargs.get('P_M_M_interaction_limit',1)
        # try_smaller_steps = kwargs.get('try_smaller_steps', True)
        # control_dir=kwargs.get('control_dir','L')  # L for lateral and V for Vertical
        # ops_analysis=kwargs.get('analysis','proportional_limit_point')
        lateral_load_scale=kwargs.get('lateral_load_scale',1)
        vertical_load_scale=kwargs.get('vertical_load_scale',1)
        plot=kwargs.get('plot')
        plot_defo=kwargs.get('plot_defo',False)
        analysis_msg=kwargs.get("analysis_msg",0)


        # Initialize analysis results
        results = AnalysisResults()
        attributes = ['load_ratio','vertical_reaction','base_shear','control_node_displacement', 'control_node_displacement_absolute',
                      'lowest_eigenvalue','absolute_maximum_strain','max_P_M_M_interaction','P_M_M_interaction_all_elements','Element_Forces']
        
        for attr in attributes:
            setattr(results, attr, [])

        # Define function to find limit point
        def find_limit_point():
            if Structures_2D.print_ops_status:
                print(results.exit_message)
            if 'Moving to Displacement Controlled Analysis' in results.exit_message:
                return
            if 'Analysis Failed' == results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif 'Eigenvalue Limit Reached' == results.exit_message:
                ind, x = find_limit_point_in_list(results.lowest_eigenvalue, eigenvalue_limit)
            elif 'Extreme Steel Fiber Strain Limit Reached' == results.exit_message:
                ind, x = find_limit_point_in_list(results.absolute_maximum_strain, steel_strain_limit)
            elif 'P_M_M interaction Limit Reached' == results.exit_message:
                ind, x = find_limit_point_in_list(results.max_P_M_M_interaction, P_M_M_interaction_limit)            
            else:
                raise Exception('Unknown limit point')
            results.maximum_load_ratio_at_limit_point = interpolate_list(results.load_ratio, ind, x)
            print('inside find limit point, Max Load Ratio',results.maximum_load_ratio_at_limit_point)
            

        def plot_analysis_history():

            line_plot(
                results.control_node_displacement,
                results.load_ratio,
                xlabel="Displacement at Control Node",
                ylabel="Load Ratio λ",
                title="Load Ratio vs Displacement",
                show=True,
            )

            line_plot(
                results.load_ratio,
                results.lowest_eigenvalue,
                xlabel="Load Ratio λ",
                ylabel="Lowest Eigenvalue",
                title="Lowest Eigenvalue vs Load Ratio",
                show=True,
            )

            line_plot(
                results.vertical_reaction,
                results.load_ratio,
                xlabel="Vertical Reaction",
                ylabel="Load Ratio λ",
                title="Load Ratio vs Vertical Reaction",
                show=True,
            )

            line_plot(
                results.base_shear,
                results.load_ratio,
                xlabel="Base Shear",
                ylabel="Load Ratio λ",
                title="Load Ratio vs Base Shear",
                show=True,
            )

            line_plot(
                results.absolute_maximum_strain,
                results.load_ratio,
                xlabel="Absolute Maximum Strain",
                ylabel="Load Ratio λ",
                title="Load Ratio vs Absolute Maximum Strain",
                show=True,
            )

            line_plot(
                results.load_ratio,
                results.max_P_M_M_interaction,
                xlabel="Load Ratio λ",
                ylabel="Max P-M-M Interaction",
                title="Max P-M-M Interaction vs Load Ratio",
                show=True,
            )
            
        fail_during_LCA=True

        ops.timeSeries('Linear', self.load_timeseries_counter)
        ops.pattern('Plain',self.load_pattern_counter, self.load_timeseries_counter)
        self.add_vertical_dead_live_wall_notional_loads(vertical_load_scale=vertical_load_scale)
        self.add_lateral_wind_loads(lateral_load_scale=lateral_load_scale)
        # region Define recorder
        def record():
            time = ops.getTime()
            results.load_ratio.append(time)
            ops.reactions()
            total_vertical_rxn=sum(ops.nodeReaction(n[0])[1] for n in self.NODES_TO_FIX)
            base_shear=sum(ops.nodeReaction(n[0])[0] for n in self.NODES_TO_FIX)
            results.vertical_reaction.append(total_vertical_rxn)
            results.base_shear.append(base_shear)
            results.lowest_eigenvalue.append(ops.eigen("-genBandArpack", 1)[0])
            results.absolute_maximum_strain.append(self.return_max_of_fiber_strain_in_all_elements())
            results.control_node_displacement.append(ops.nodeDisp(control_node, control_dof))
            max_PMM, max_ele_tag,P_M_M_interaction_all_elements,Element_Forces=self.return_P_M_M_interaction_values()
            results.max_P_M_M_interaction.append(max_PMM)
            results.P_M_M_interaction_all_elements.append(P_M_M_interaction_all_elements)
            results.Element_Forces.append(Element_Forces)
        # endregion
        control_node,control_dof=self.get_control_node_and_dof(control_dir='L')
        # Create output folder
        os.makedirs(os.path.join("Column_Results/Analysis_History", self.Frame_id), exist_ok=True)

        ops.constraints('Plain')
        # ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-8, 10,analysis_msg)
        ops.algorithm('Newton')
        ops.integrator('LoadControl',1/num_steps_LCA) 
        ops.analysis('Static')

        record()
        for i in range(num_steps_LCA):
            if Structures_2D.print_ops_status:
                print(f'Running Load Controlled Analysis Step {i}')
            ok = ops.analyze(1)
            if ok != 0:
                print(f'Load controlled analysis failed in step {i}')
                results.exit_message = 'Analysis Failed'
                find_limit_point()
                drift=self.return_drift_of_all_storeys_at_given_axis(1)
                if plot:
                    plot_analysis_history() 
                return drift,results,fail_during_LCA
            else:
                print('Load controlled analysis PASSED')
            record()

            # Check for lowest eigenvalue less than zero
            # if eigenvalue_limit is not None:
            #     if results.lowest_eigenvalue[-1] < eigenvalue_limit:
            #         results.exit_message = 'Eigenvalue Limit Reached'
            #         drift=self.return_drift_of_all_storeys_at_given_axis(1)
            #         find_limit_point()
            #         if plot:
            #             plot_analysis_history()
            #         return drift,results,fail_during_LCA
                    # break

            # Check for strain in extreme steel fiber
            # if steel_strain_limit is not None:
            #     # if Structures_2D.print_ops_status:
            #     #     print(f'Checking Steel Tensile Strain')
            #     if results.absolute_maximum_strain[-1] > steel_strain_limit:
            #         results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
            #         drift=self.return_drift_of_all_storeys_at_given_axis(1)
            #         find_limit_point()
            #         if plot:
            #             plot_analysis_history()
            #         return drift,results,fail_during_LCA
                    # break
            # Check for maximum PMM interaction value    

            # if self.Elastic_analysis:
            #     if P_M_M_interaction_limit is not None:
            #         # if Structures_2D.print_ops_status:
            #         #     print(f'Checking PMM Interaction')
            #         if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
            #             results.exit_message = 'P_M_M interaction Limit Reached'
            #             drift=self.return_drift_of_all_storeys_at_given_axis(1)
            #             find_limit_point()
            #             if plot:
            #                 plot_analysis_history()
            #             return drift,results,fail_during_LCA

        ## Store drift at all storeys at axis 1 ##
        
        drift=self.return_drift_of_all_storeys_at_given_axis(1)
        print(results.control_node_displacement)
        print(results.lowest_eigenvalue)
        if plot:
            plot_analysis_history()
            plt.show()
            

        if plot_defo:
            try:
                # import opsvis
                opsvis.plot_defo()
                plt.show()
                
            except:
                print("opsvis not available for deformation plotting.")


        return drift,results,fail_during_LCA

    def get_lateral_loading_direction(self):
        print('I am inside get_lateral_loading_direction method')
        # input()
        # if self.geometric_imperfection_ratio>0:
        #     self.wind_load_dirn='right'
        # elif self.geometric_imperfection_ratio<0:
        #     self.wind_load_dirn='left'
        # else:
        if self.wind_load_dirn is None:
            Dummy_Frame=self.rebuild_with_overrides(Second_order_effects=True,
                                                        Residual_Stress=False,
                                                        Elastic_analysis=True,
                                                        stiffness_reduction=0.8,
                                                        strength_reduction=1,
                                                        Geometric_Imperfection=False,
                                                        geometric_imperfection_ratio=1/500,
                                                        Notional_load=False)
            Dummy_Frame.generate_Nodes_and_Element_Connectivity()
            Dummy_Frame.build_ops_model()
            _,dummy_results,_=Dummy_Frame.run_load_controlled_anlaysis(lateral_load_scale=0.001,vertical_load_scale=0.1,plot=False)
            # if dummy_results.lowest_eigenvalue[-1]<0:
            #     displacement_to_check=dummy_results.control_node_displacement[-2]
            # else:
            displacement_to_check=dummy_results.control_node_displacement[-1]

            print(dummy_results.control_node_displacement)
            print(dummy_results.lowest_eigenvalue)
            
            if displacement_to_check>=0:
                self.wind_load_dirn='right'
            else:
                self.wind_load_dirn='left'
            ops.wipe()

        return self.wind_load_dirn  

    def get_initial_out_of_straightness_direction(self):
        print('I am inside get_initial_out_of_straightness_direction method')

        if self.initial_out_of_straightness_dirn is None:
            print("Right out of straightness model running")
            Dummy_Frame_right=self.rebuild_with_overrides(Second_order_effects=True,
                                                        Residual_Stress=True,
                                                        Elastic_analysis=False,
                                                        stiffness_reduction=0.9,
                                                        strength_reduction=0.9,
                                                        Geometric_Imperfection=True,
                                                        initial_out_of_straightness_dirn='right',
                                                        wind_load_dirn=self.wind_load_dirn,
                                                        Notional_load=False)
            Dummy_Frame_right.generate_Nodes_and_Element_Connectivity()
            Dummy_Frame_right.create_distorted_nodes_and_element_connectivity()
            Dummy_Frame_right.build_ops_model()



            # opsv.plot_model()
            # opsv.plot_load()
            # Dummy_Frame_right.plot_model()
            

            dummy_results_right,fail_during_LCA =Dummy_Frame_right.run_displacement_controlled_analysis(target_disp=1,steps=1000,plot_defo=False,num_steps_LCA=50, 
                                                      analysis='proportional_limit_point',
                                                      vertical_load_scale=1,
                                                      lateral_load_scale=0.5,
                                                      control_dir='L',try_smaller_steps=False,
                                                      live_plot=False,tolerance=1e-6,iterations=10) 



            print("Left out of straightness model running")
            Dummy_Frame_left=self.rebuild_with_overrides(Second_order_effects=True,
                                                        Residual_Stress=True,
                                                        Elastic_analysis=False,
                                                        stiffness_reduction=0.9,
                                                        strength_reduction=0.9,
                                                        Geometric_Imperfection=True,
                                                        initial_out_of_straightness_dirn='left',
                                                        wind_load_dirn=self.wind_load_dirn,
                                                        Notional_load=False)
            Dummy_Frame_left.generate_Nodes_and_Element_Connectivity()
            Dummy_Frame_left.create_distorted_nodes_and_element_connectivity()
            Dummy_Frame_left.build_ops_model()

            # opsv.plot_model()
            # opsv.plot_load()
            # Dummy_Frame_left.plot_model()

            
            dummy_results_left,fail_during_LCA =Dummy_Frame_left.run_displacement_controlled_analysis(target_disp=1,steps=1000,plot_defo=False,num_steps_LCA=50, 
                                                      analysis='proportional_limit_point',
                                                      vertical_load_scale=1,
                                                      lateral_load_scale=0.5,
                                                      control_dir='L',try_smaller_steps=False,
                                                      live_plot=False,tolerance=1e-6,iterations=10) 


            print(dummy_results_right.maximum_load_ratio_at_limit_point)
            print(dummy_results_left.maximum_load_ratio_at_limit_point)
            if dummy_results_right.maximum_load_ratio_at_limit_point>dummy_results_left.maximum_load_ratio_at_limit_point:
                self.initial_out_of_straightness_dirn='left'
            else:
                self.initial_out_of_straightness_dirn='right'
            ops.wipe()

        return self.initial_out_of_straightness_dirn    
                 
    def get_del2_over_del1(self,vertical_load_scale=1, lateral_load_scale=1):
        print('Second order Frame')
        # print(self)
        Second_order_frame=self.rebuild_with_overrides(Second_order_effects=True,
                                                       Residual_Stress=False,
                                                       Elastic_analysis=True,
                                                       stiffness_reduction=0.8,
                                                       strength_reduction=1,
                                                       Geometric_Imperfection=False,
                                                       geometric_imperfection_ratio=1/500,
                                                       Notional_load=False,
                                                       wind_load_dirn=self.wind_load_dirn)
        print(Second_order_frame.Second_order_effects)
        print(Second_order_frame.Geometric_Imperfection)
        print(Second_order_frame.geometric_imperfection_ratio)
        print(Second_order_frame.Notional_load)
        print(Second_order_frame.stiffness_reduction)
        print(Second_order_frame.strength_reduction)
        print(Second_order_frame.Residual_Stress)
        print(Second_order_frame.Elastic_analysis)
        
        Second_order_frame.generate_Nodes_and_Element_Connectivity()
        Second_order_frame.create_distorted_nodes_and_element_connectivity()
        Second_order_frame.build_ops_model()
        drift_with_second_order_effects,results,fail_during_LCA=Second_order_frame.run_load_controlled_anlaysis(vertical_load_scale=vertical_load_scale,lateral_load_scale=lateral_load_scale,plot=False)

        if results.lowest_eigenvalue[-1]<0:
            print('Warning: The lowest eigenvalue is negative, which indicates that the structure has gone past elastic critical buckling limit. So the drift is treated as infinite.')
            drift_with_second_order_effects = [float('inf')] * len(drift_with_second_order_effects)
        print('drift_with_second_order_effects',drift_with_second_order_effects)
        # print(self)

        # input('Press Enter to continue...')
        print('First order Frame')
        First_order_frame=self.rebuild_with_overrides(Second_order_effects=False,
                                                       Residual_Stress=False,
                                                       Elastic_analysis=True,
                                                       stiffness_reduction=0.8,
                                                       strength_reduction=1,
                                                       Geometric_Imperfection=False,
                                                       geometric_imperfection_ratio=1/500,
                                                       Notional_load=False,
                                                       wind_load_dirn=self.wind_load_dirn)
        print(First_order_frame.Second_order_effects)
        print(First_order_frame.Geometric_Imperfection)
        print(First_order_frame.geometric_imperfection_ratio)
        print(First_order_frame.Notional_load)
        print(First_order_frame.stiffness_reduction)
        print(First_order_frame.strength_reduction)
        print(First_order_frame.Residual_Stress)
        print(First_order_frame.Elastic_analysis)
        
        First_order_frame.generate_Nodes_and_Element_Connectivity()
        First_order_frame.create_distorted_nodes_and_element_connectivity()
        First_order_frame.build_ops_model()
        drift_with_first_order_effects,results,fail_during_LCA=First_order_frame.run_load_controlled_anlaysis(vertical_load_scale=vertical_load_scale,lateral_load_scale=lateral_load_scale,plot=False)
        if results.lowest_eigenvalue[-1]<0:
            raise Exception('Warning: Got negative eigenvalue for 1st order elastic analysis while calculating del2_over_del1.')        
        print('drift_with_first_order_effects',drift_with_first_order_effects)
        # input('Press Enter to continue...')

        del2_over_del1=[]
        for i in range(len(drift_with_first_order_effects)):
            del2_over_del1.append(drift_with_second_order_effects[i]/drift_with_first_order_effects[i])
        print('del2_over_del1',del2_over_del1)
        
        return max(del2_over_del1)


    def run_displacement_controlled_analysis(self, target_disp=1, steps=10000, plot_defo=False,**kwargs):
        """
        Runs displacement-controlled analysis and plots load ratio (λ) vs. displacement and vertical reaction.

        Parameters:
            target_disp (float): Target horizontal displacement at control node 
            steps (int): Number of steps to reach target
            plot_defo (bool): Whether to plot deformed shape at end
        """
        
        incr_LCA= kwargs.get('incr_LCA', 0.01)          ######### LCA refers to Load Controlled Analysis
        num_steps_LCA= kwargs.get('num_steps_LCA', 50)            ######### LCA refers to Load Controlled Analysis
        steel_strain_limit = kwargs.get('steel_strain_limit', 0.05)
        eigenvalue_limit = kwargs.get('eigenvalue_limit', 0)
        P_M_M_interaction_limit=kwargs.get('P_M_M_interaction_limit',1)
        try_smaller_steps = kwargs.get('try_smaller_steps', False)
        control_dir=kwargs.get('control_dir','L')  # L for lateral and V for Vertical
        ops_analysis=kwargs.get('analysis','proportional_limit_point')
        lateral_load_scale=kwargs.get('lateral_load_scale',1)
        vertical_load_scale=kwargs.get('vertical_load_scale',1)
        live_plot = kwargs.get("live_plot", False)
        live_plot_every = max(1, int(kwargs.get("live_plot_every", 1)))
        tol=kwargs.get("tolerance", 1e-10)
        iter=kwargs.get("iterations", 10)
        analysis_msg=kwargs.get("analysis_msg",0)


        # Initialize analysis results
        results = AnalysisResults()
        attributes = ['load_ratio','vertical_reaction','base_shear','control_node_displacement', 'control_node_displacement_absolute',
                      'lowest_eigenvalue','absolute_maximum_strain','max_P_M_M_interaction','P_M_M_interaction_all_elements','Element_Forces']
        
        def initialize_results():
            for attr in attributes:
                setattr(results, attr, [])

        initialize_results()

        live_fig = None
        live_axes = None
        live_lines = {}
        was_interactive = plt.isinteractive()

        if live_plot:
            plt.ion()

            live_fig, live_axes = plt.subplots(
                nrows=2,
                ncols=3,
                figsize=(16, 9),
                facecolor="white",
            )

            axes = live_axes.ravel()

            # Store both the history lines and current-point markers
            live_lines = {}
            live_last_points = {}

            subplot_settings = [
                {
                    "key": "displacement",
                    "xlabel": "Displacement at Control Node",
                    "ylabel": "Load Ratio λ",
                    "title": "Load Ratio vs Displacement",
                },
                {
                    "key": "eigenvalue",
                    "xlabel": "Load Ratio λ",
                    "ylabel": "Lowest Eigenvalue",
                    "title": "Lowest Eigenvalue vs Load Ratio",
                },
                {
                    "key": "vertical_reaction",
                    "xlabel": "Vertical Reaction",
                    "ylabel": "Load Ratio λ",
                    "title": "Load Ratio vs Vertical Reaction",
                },
                {
                    "key": "base_shear",
                    "xlabel": "Base Shear",
                    "ylabel": "Load Ratio λ",
                    "title": "Load Ratio vs Base Shear",
                },
                {
                    "key": "strain",
                    "xlabel": "Absolute Maximum Strain",
                    "ylabel": "Load Ratio λ",
                    "title": "Load Ratio vs Maximum Strain",
                },
                {
                    "key": "pmm",
                    "xlabel": "Load Ratio λ",
                    "ylabel": "Maximum P-M-M Interaction",
                    "title": "P-M-M Interaction vs Load Ratio",
                },
            ]

            for ax, settings in zip(axes, subplot_settings):
                key = settings["key"]

                # Complete analysis history
                live_lines[key], = ax.plot(
                    [],
                    [],
                    color="tab:blue",
                    linewidth=2.0,
                    marker="o",
                    markersize=6,
                    markevery=1,
                    alpha=0.9,
                    zorder=2,
                )

                # Current/last point
                live_last_points[key], = ax.plot(
                    [],
                    [],
                    marker="o",
                    linestyle="None",
                    markersize=6,
                    markerfacecolor="red",
                    markeredgecolor="black",
                    markeredgewidth=0.8,
                    zorder=5,
                    label="Current point",
                )

                ax.set_xlabel(
                    settings["xlabel"],
                    fontsize=10,
                    fontweight="medium",
                )
                ax.set_ylabel(
                    settings["ylabel"],
                    fontsize=10,
                    fontweight="medium",
                )
                ax.set_title(
                    settings["title"],
                    fontsize=11,
                    fontweight="bold",
                    pad=10,
                )

                ax.grid(
                    True,
                    linestyle="--",
                    linewidth=0.7,
                    alpha=0.4,
                )

                ax.tick_params(
                    axis="both",
                    labelsize=9,
                    direction="in",
                )

                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)

                ax.margins(x=0.05, y=0.08)

            # Only one legend is needed
            axes[0].legend(
                loc="best",
                fontsize=9,
                frameon=True,
            )

            live_fig.suptitle(
                f"Live Analysis History: {self.Frame_id}",
                fontsize=16,
                fontweight="bold",
                y=0.98,
            )

            live_fig.tight_layout(
                rect=[0.02, 0.02, 0.98, 0.95],
                h_pad=2.0,
                w_pad=2.0,
            )

            plt.show(block=False)
            plt.pause(0.1)

        def update_live_plots(force=False):
            if not live_plot:
                return

            if live_fig is None:
                return

            if not plt.fignum_exists(live_fig.number):
                return

            number_of_points = len(results.load_ratio)

            if number_of_points == 0:
                return

            if not force and number_of_points % live_plot_every != 0:
                return

            plot_data = {
                "displacement": (
                    results.control_node_displacement,
                    results.load_ratio,
                ),
                "eigenvalue": (
                    results.load_ratio,
                    results.lowest_eigenvalue,
                ),
                "vertical_reaction": (
                    results.vertical_reaction,
                    results.load_ratio,
                ),
                "base_shear": (
                    results.base_shear,
                    results.load_ratio,
                ),
                "strain": (
                    results.absolute_maximum_strain,
                    results.load_ratio,
                ),
                "pmm": (
                    results.load_ratio,
                    results.max_P_M_M_interaction,
                ),
            }

            for key, (x_values, y_values) in plot_data.items():
                if len(x_values) == 0 or len(y_values) == 0:
                    continue

                # Update the complete curve
                live_lines[key].set_data(
                    x_values,
                    y_values,
                )

                # Highlight the newest point
                live_last_points[key].set_data(
                    [x_values[-1]],
                    [y_values[-1]],
                )

            for ax in live_axes.ravel():
                ax.relim()
                ax.autoscale_view()

            live_fig.canvas.draw_idle()
            live_fig.canvas.flush_events()

            plt.pause(0.005)
        
        
        # Define function to find limit point
        def find_limit_point():

            
            print(results.exit_message)

            if 'Moving to Displacement Controlled Analysis' == results.exit_message:
                return
            elif 'Analysis Failed in Displacement Controlled Loading in Non Proportional Analysis'==results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif  'Analysis Failed In Load Controlled Loading before entering Displacement controlled Loading' == results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif 'Analysis Failed' == results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif 'Eigenvalue Limit Reached' == results.exit_message:
                ind, x = find_limit_point_in_list(results.lowest_eigenvalue, eigenvalue_limit)
            elif 'Extreme Steel Fiber Strain Limit Reached' == results.exit_message:
                ind, x = find_limit_point_in_list(results.absolute_maximum_strain, steel_strain_limit)
            elif 'P_M_M interaction Limit Reached' == results.exit_message:
                ind, x = find_limit_point_in_list(results.max_P_M_M_interaction, P_M_M_interaction_limit)            
            else:
                raise Exception('Unknown limit point')
            results.maximum_load_ratio_at_limit_point = interpolate_list(results.load_ratio, ind, x)
            print('Line 884, Max Load Ratio',results.maximum_load_ratio_at_limit_point)


        if ops_analysis.lower()=='proportional_limit_point':

            fail_during_LCA=True

            ops.timeSeries('Linear', self.load_timeseries_counter)
            ops.pattern('Plain',self.load_pattern_counter, self.load_timeseries_counter)
            self.add_vertical_dead_live_wall_notional_loads(vertical_load_scale=vertical_load_scale)
            self.add_lateral_wind_loads(lateral_load_scale=lateral_load_scale)
            # region Define recorder
            def record():
                time = ops.getTime()
                results.load_ratio.append(time)
                ops.reactions()
                total_vertical_rxn=sum(ops.nodeReaction(n[0])[1] for n in self.NODES_TO_FIX)
                base_shear=sum(ops.nodeReaction(n[0])[0] for n in self.NODES_TO_FIX)
                results.vertical_reaction.append(total_vertical_rxn)
                results.base_shear.append(base_shear)
                results.lowest_eigenvalue.append(ops.eigen("-genBandArpack", 1)[0])
                results.absolute_maximum_strain.append(self.return_max_of_fiber_strain_in_all_elements())
                results.control_node_displacement.append(ops.nodeDisp(control_node, control_dof))
                max_PMM, max_ele_tag,P_M_M_interaction_all_elements,Element_Forces=self.return_P_M_M_interaction_values()
                results.max_P_M_M_interaction.append(max_PMM)
                results.P_M_M_interaction_all_elements.append(P_M_M_interaction_all_elements)
                results.Element_Forces.append(Element_Forces)
            # endregion

            control_node,control_dof=self.get_control_node_and_dof(control_dir=control_dir)

            # Create output folder
            os.makedirs(os.path.join("Column_Results/Analysis_History", self.Frame_id), exist_ok=True)

            ops.constraints('Plain')
            # ops.constraints('Transformation')
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test('NormUnbalance', tol, iter, analysis_msg)
            ops.algorithm('RaphsonNewton')  
            ops.analysis('Static')
            dU = target_disp / steps

            print(self.wind_load_dirn)
            
            if self.wind_load_dirn=="left" or control_dof==2:
                dU=-dU 
            ops.integrator('DisplacementControl', control_node, control_dof, dU)

            record()

            update_live_plots()

            i=1
            # save_deformed_shape_frame(i,"temporary_folder")
            while True:
                print(f'Running Displacement Controlled Analysis {i}')
                i=i+1
                fail_during_LCA=False
                ok = ops.analyze(1)

                if try_smaller_steps:
                    if ok != 0:
                        if Structures_2D.print_ops_status:
                            print(f'Failed for Step size of {dU}. Now Trying the step size of: {dU / 10}')
                        ops.integrator('DisplacementControl',control_node, control_dof, dU / 10)
                        ok = ops.analyze(1)

                    if ok != 0:
                        if Structures_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 100}')
                        ops.integrator('DisplacementControl',control_node, control_dof, dU / 100)
                        ok = ops.analyze(1)

                    if ok != 0:
                        if Structures_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 1000}')
                        ops.integrator('DisplacementControl', control_node, control_dof, dU / 1000)
                        ok = ops.analyze(1)
                        if ok == 0:
                            # dU = dU / 10
                            if Structures_2D.print_ops_status:
                                print(f'Changed the step size to: {dU}')

                    if ok != 0:
                        if Structures_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 10000}')
                        ops.integrator('DisplacementControl', control_node, control_dof, dU / 10000)
                        ok = ops.analyze(1)
                        if ok == 0:
                            # dU = dU / 10
                            if Structures_2D.print_ops_status:
                                print(f'Changed the step size to: {dU }')

                if ok != 0:
                    if Structures_2D.print_ops_status:
                        print('Trying ModifiedNewton')
                    ops.algorithm('ModifiedNewton')
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Structures_2D.print_ops_status:
                            print('ModifiedNewton worked')

                if ok != 0:
                    if Structures_2D.print_ops_status:
                        print('Trying KrylovNewton')
                    ops.algorithm('KrylovNewton')
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Structures_2D.print_ops_status:
                            print('KrylovNewton worked')

                if ok != 0:
                    if Structures_2D.print_ops_status:
                        print('Trying KrylovNewton and Greater Tolerance')
                    ops.algorithm('KrylovNewton')
                    ops.test('NormUnbalance', tol*100, iter,analysis_msg)
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Structures_2D.print_ops_status:
                            print('KrylovNewton worked')

                if ok == 0:
                    # Reset analysis options
                    print(f'Displacement controlled analysis step {i-1} PASSED')
                    ops.algorithm('RaphsonNewton')
                    ops.test('NormUnbalance', tol, iter,analysis_msg)
                    ops.integrator('DisplacementControl', control_node, control_dof, dU)
                else:
                    print('Analysis Failed')
                    results.exit_message = 'Analysis Failed'
                    break


                record()
                update_live_plots()


                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    print(results.lowest_eigenvalue[-1])
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        break

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    # if Structures_2D.print_ops_status:
                    #     print(f'Checking Steel Tensile Strain')
                    print(results.absolute_maximum_strain[-1])
                    if results.absolute_maximum_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                        break
                # Check for maximum PMM interaction value    
                if self.Elastic_analysis:
                    if P_M_M_interaction_limit is not None:
                        # if Structures_2D.print_ops_status:
                        #     print(f'Checking PMM Interaction')
                        if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                            results.exit_message = 'P_M_M interaction Limit Reached'
                            break
                
                

                    
                # if i%(steps/20)==0:
                     # save_deformed_shape_frame(i,"temporary_folder")
            update_live_plots(force=True)

            if live_plot and not was_interactive:
                plt.ioff()

            # save_deformed_shape_frame(i,"temporary_folder")

            find_limit_point()


        elif ops_analysis.lower()=='non_proportional_limit_point':
            
            fail_during_LCA=True

            ops.timeSeries('Linear', self.load_timeseries_counter)
            ops.pattern('Plain',self.load_pattern_counter, self.load_timeseries_counter)
            self.add_vertical_dead_live_wall_notional_loads(vertical_load_scale=vertical_load_scale)

            # region Define recorder
            def record():
                time = ops.getTime()
                results.load_ratio.append(time)
                ops.reactions()
                total_vertical_rxn=sum(ops.nodeReaction(n[0])[1] for n in self.NODES_TO_FIX)
                base_shear=sum(ops.nodeReaction(n[0])[0] for n in self.NODES_TO_FIX)
                results.vertical_reaction.append(total_vertical_rxn)
                results.base_shear.append(base_shear)
                results.lowest_eigenvalue.append(ops.eigen("-genBandArpack", 1)[0])
                results.absolute_maximum_strain.append(self.return_max_of_fiber_strain_in_all_elements())
                results.control_node_displacement.append(ops.nodeDisp(control_node, control_dof))
                max_PMM, max_ele_tag,P_M_M_interaction_all_elements,Element_Forces=self.return_P_M_M_interaction_values()
                results.max_P_M_M_interaction.append(max_PMM)
                results.P_M_M_interaction_all_elements.append(P_M_M_interaction_all_elements)
                results.Element_Forces.append(Element_Forces)
            # endregion

            control_node,control_dof=self.get_control_node_and_dof(control_dir=control_dir)


            # Create output folder
            os.makedirs(os.path.join("Column_Results/Analysis_History", self.Frame_id), exist_ok=True)

            # ops.constraints('Plain')
            ops.constraints('Plain')
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test('NormUnbalance', 1e-8, iter,analysis_msg)
            ops.algorithm('Newton')
            ops.integrator('LoadControl',1/num_steps_LCA)  ## 1/num_steps_LCA because we want the entire vertical load to be applied before starting displacement controlled analysis.
            ops.analysis('Static')
            record()
            for i in range(num_steps_LCA):
                if Structures_2D.print_ops_status:
                    print(f'Running Load Controlled Analysis Step {i}')
                ok = ops.analyze(1)

                if ok != 0:
                    print(f'Load controlled analysis failed in step {i}')
                    results.exit_message = 'Analysis Failed In Load Controlled Loading before entering Displacement controlled Loading'
                    find_limit_point()
                    return results,fail_during_LCA
                else:
                    print('Load controlled analysis PASSED')
                    results.exit_message='Moving to Displacement Controlled Analysis'
                record()
                update_live_plots()


                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        find_limit_point()
                        return results,fail_during_LCA
                        # break

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    # if Structures_2D.print_ops_status:
                    #     print(f'Checking Steel Tensile Strain')
                    if results.absolute_maximum_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                        find_limit_point()
                        return results,fail_during_LCA
                        # break
                # Check for maximum PMM interaction value    
                if self.Elastic_analysis:
                    if P_M_M_interaction_limit is not None:
                        # if Structures_2D.print_ops_status:
                        #     print(f'Checking PMM Interaction')
                        if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                            results.exit_message = 'P_M_M interaction Limit Reached'
                            find_limit_point()
                            return results,fail_during_LCA
                            # break
            # find_limit_point()


            dU = target_disp / steps
            if self.wind_load_dirn=="left" or control_dof==2:
                dU=-dU 

            ops.loadConst('-time', 0.0)
            ops.timeSeries('Linear', self.load_timeseries_counter+1)
            ops.pattern('Plain',self.load_pattern_counter+1, self.load_timeseries_counter+1)
            self.add_lateral_wind_loads(lateral_load_scale=lateral_load_scale)
            ops.algorithm('RaphsonNewton')
            ops.test('NormUnbalance', tol, iter,analysis_msg)
            ops.integrator('DisplacementControl', control_node, control_dof, dU)

            record()

            update_live_plots()

            initialize_results()

            i=1
            while True:
                print(f'Running Displacement Controlled Analysis {i}')
                i=i+1
                fail_during_LCA=False
                ok = ops.analyze(1)
                if try_smaller_steps:
                    if ok != 0:
                        if Structures_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 10}')
                        ops.integrator('DisplacementControl',control_node, control_dof, dU / 10)
                        ok = ops.analyze(1)

                    if ok != 0:
                        if Structures_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 100}')
                        ops.integrator('DisplacementControl',control_node, control_dof, dU / 100)
                        ok = ops.analyze(1)

                    if ok != 0:
                        if Structures_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 1000}')
                        ops.integrator('DisplacementControl', control_node, control_dof, dU / 1000)
                        ok = ops.analyze(1)
                        if ok == 0:
                            # dU = dU / 10
                            if Structures_2D.print_ops_status:
                                print(f'Changed the step size to: {dU}')

                    if ok != 0:
                        if Structures_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 10000}')
                        ops.integrator('DisplacementControl', control_node, control_dof, dU / 10000)
                        ok = ops.analyze(1)
                        if ok == 0:
                            # dU = dU / 10
                            if Structures_2D.print_ops_status:
                                print(f'Changed the step size to: {dU / 10}')

                if ok != 0:
                    if Structures_2D.print_ops_status:
                        print('Trying ModifiedNewton')
                    ops.algorithm('ModifiedNewton')
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Structures_2D.print_ops_status:
                            print('ModifiedNewton worked')

                if ok != 0:
                    if Structures_2D.print_ops_status:
                        print('Trying KrylovNewton')
                    ops.algorithm('KrylovNewton')
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Structures_2D.print_ops_status:
                            print('KrylovNewton worked')

                if ok != 0:
                    if Structures_2D.print_ops_status:
                        print('Trying KrylovNewton and Greater Tolerance')
                    ops.algorithm('KrylovNewton')
                    ops.test('NormUnbalance', tol*100, iter,analysis_msg)
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Structures_2D.print_ops_status:
                            print('KrylovNewton worked')

                if ok == 0:
                    # Reset analysis options
                    print(f'Displacement controlled analysis step {i-1} PASSED')
                    ops.algorithm('RaphsonNewton')
                    ops.test('NormUnbalance', tol, iter,analysis_msg)
                    ops.integrator('DisplacementControl', control_node, control_dof, dU)
                else:
                    print('Analysis Failed in Displacement Controlled Loading in Non Proportional Analysis')
                    results.exit_message = 'Analysis Failed in Displacement Controlled Loading in Non Proportional Analysis'
                    break


                record()
                update_live_plots()


                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        break

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    # if Structures_2D.print_ops_status:
                    #     print(f'Checking Steel Tensile Strain')
                    if results.absolute_maximum_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                        break
                # Check for maximum PMM interaction value    
                if self.Elastic_analysis:
                    if P_M_M_interaction_limit is not None:
                        # if Structures_2D.print_ops_status:
                        #     print(f'Checking PMM Interaction')
                        if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                            results.exit_message = 'P_M_M interaction Limit Reached'
                            break

            update_live_plots(force=True)
            find_limit_point()
       
            if live_plot and not was_interactive:
                plt.ioff()

        else:
            raise Exception('Give valid ops_analysis option')
        # Optional: plot deformed shape
        if plot_defo:
            try:
                # import opsvis
                opsvis.plot_defo()
            except:
                print("opsvis not available for deformation plotting.")

        return results,fail_during_LCA



    def save_moments_by_member(self, filename='max_member_moments.csv'):
        os.makedirs(self.Frame_id, exist_ok=True)
        full_path = os.path.join(self.Frame_id, filename)

        conversion_factor = 8.85074579  # kN·m to kip·in
        member_moment_data = []

        for member in self.member_list:
            member_tag = member[0]
            element_tags = member[1:]

            max_moment = 0.0

            for eleTag in element_tags:
                try:
                    forces = ops.eleResponse(eleTag, 'localForce')
                    Mz_i_kNm = forces[2]
                    Mz_j_kNm = forces[5]

                    # Convert to kip-in
                    Mz_i = Mz_i_kNm * conversion_factor
                    Mz_j = Mz_j_kNm * conversion_factor

                    max_moment = max(max_moment, abs(Mz_i), abs(Mz_j))

                except Exception as e:
                    print(f"Error in element {eleTag} of member {member_tag}: {e}")

            member_moment_data.append({
                'Member': member_tag,
                'Max_Abs_Moment (kip-in)': max_moment
            })

        df = pd.DataFrame(member_moment_data)
        df.to_csv(full_path, index=False)
        print(f"Saved member-level max moments (kip-in) to: {full_path}")

    def plot_all_fiber_section_in_the_model(self):
        for sec_tag in self.beam_section_tags.values():
            get_fiber_data(f'{sec_tag}',plot_fibers=True)

        for sec_tag,_ in self.column_section_tags.values():
            get_fiber_data(f'{sec_tag}',plot_fibers=True)

    def reset_analysis():
        """
        Resets the analysis by setting time to 0,
        removing the recorders and wiping the analysis.
        """    
        
        # Reset for next analysis case
        ##  Set the time in the Domain to zero
        ops.setTime(0.0)
        ## Set the loads constant in the domain
        ops.loadConst()
        ## Remove all recorder objects.
        ops.remove('recorders')
        ## destroy all components of the Analysis object
        ops.wipeAnalysis()
        
    
    def plot_model(self):
        plot_undeformed_2d(axis_equal=True)
        

    def display_node_coords(self):
        get_node_coords_and_disp()

    def plot_deformed_shape(self):
        plot_deformed_2d(axis_equal=True,scale_factor=8)


    