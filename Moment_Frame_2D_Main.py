import numpy as np
import math as math
import pandas as pd
import openseespy.opensees as ops
import os
import opsvis as opsv
from libdenavit.OpenSees.plotting import *
from libdenavit.OpenSees.get_fiber_data import *
from Units import *
from math import pi, ceil
from libdenavit.section.wide_flange import *
from libdenavit.OpenSees import AnalysisResults
from libdenavit import find_limit_point_in_list, interpolate_list
import copy

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

class Moment_Frame_2D:

    print_ops_status=True

    def __init__(self,width_of_bay,storey_height,
                no_of_elements_column, no_of_elements_beam,
                 beam_section,column_section,load_combination_multipliers,Frame_id,
                 **kwargs):

        self.bay_width=[0]+width_of_bay
        self.storey_height = [0] + storey_height
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
                  'Second_order_effects':True,
                  'Notional_load':False,
                  'Geometric_Imperfection':False,
                  'Residual_Stress':True,
                  'stiffness_reduction':1,
                  'strength_reduction':1,
                  'geometric_imperfection_ratio':1/500,
                  'wind_load_dirn':'right',
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
        bay = Moment_Frame_2D.nth_digit(first_node, 1)
        storey = Moment_Frame_2D.nth_digit(first_node, 2,2)
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
        axis = Moment_Frame_2D.nth_digit(first_node, 1)
        storey = Moment_Frame_2D.nth_digit(first_node, 2,2) + 1  
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

###################### The code below is written so as to change the order of elements for direct comparison with Ziemian results
        self.sorted_column_connectivity= sorted(
                        self.column_connectivity,key=lambda col: tuple(reversed(Moment_Frame_2D.find_axis_and_storey_for_columns(col[1:])))
                        )
        
        self.column_member_list = []
        group_size=self.no_of_elements_column
        for i in range(0, len(self.sorted_column_connectivity), group_size):
            group = self.sorted_column_connectivity[i:i+group_size]
            if len(group) == group_size:
                member_tag = (i // group_size) + 1
                member_element_tags = [e[0] for e in group]
                self.column_member_list.append([member_tag] + member_element_tags)


        self.sorted_beam_connectivity=sorted(self.beam_connectivity,key=lambda b: Moment_Frame_2D.find_bay_and_storey_for_beams(b[1:])[1])

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
        '''This function returns the node tags of all the floor nodes that lie in ith bay  '''
        return [
            node[0]
            for node in self.beam_nodes
            if len(str(node[0])) == 5
            and Moment_Frame_2D.nth_digit(node[0], 1) == i
            and Moment_Frame_2D.nth_digit(node[0], 2,2) != self.no_of_stories
        ]
    
    def bay_i_internal_roof_nodes(self, i):
        '''This function returns the node tags of all the roof nodes that lie in ith bay  '''
        return [
            node[0]
            for node in self.beam_nodes
            if len(str(node[0])) == 5
            and Moment_Frame_2D.nth_digit(node[0], 1) == i
            and Moment_Frame_2D.nth_digit(node[0], 2,2) == self.no_of_stories
        ]

    def axis_i_roof_nodes(self, i):
        '''This function returns the node tags of all the roof nodes that lie in ith axis  '''
        return [
            node[0]
            for node in self.beam_nodes
            if len(str(node[0])) == 3
            and Moment_Frame_2D.nth_digit(node[0], 1) == i
            and Moment_Frame_2D.nth_digit(node[0], 2,2) == self.no_of_stories
        ]
    
    def axis_i_floor_nodes(self, i):
        '''This function returns the node tags of all the floor nodes that lie in ith axis  '''
        return [
            node[0]
            for node in self.beam_nodes
            if len(str(node[0])) == 3
            and Moment_Frame_2D.nth_digit(node[0], 1) == i
            and Moment_Frame_2D.nth_digit(node[0], 2,2) != self.no_of_stories
        ]

    def bay_i_floor_load(self, i):
        ''' This function returns the point loads to be applied to all the internal floor nodes of ith bay'''
        return [((-self.D_multiplier*self.D_floor_intensity)+(-self.L_multiplier*self.L_floor_intensity))*(self.bay_width[i]/self.no_of_elements_beam)]
    
    def bay_i_roof_load(self, i):
        return [((-self.D_multiplier*self.D_roof_intensity)+(-self.L_r_multiplier*self.L_roof_intensity))*(self.bay_width[i]/self.no_of_elements_beam)]
        
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

    def create_distorted_nodes_and_element_connectivity(self,geometric_imperfection_ratio=None):
        if self.Geometric_Imperfection:
            print('Working in imperfect geometry')
            ratio = geometric_imperfection_ratio if geometric_imperfection_ratio is not None else self.geometric_imperfection_ratio
            for i in range(len(self.all_nodes)):
                    self.all_nodes[i][1] = self.all_nodes[i][1] +ratio * self.all_nodes[i][2]
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


        Steel=Steel_Material(1,E=E,Fy=Fy,G=G,Hk=Hk,density=density_of_steel)

        col_and_beam_TransTag = 1

        if self.Elastic_analysis:
            mat_type = 'Elastic'
            frc = 0
        else:
            mat_type = self.mat_type
            frc = -0.3 * Steel.Fy if self.Residual_Stress else 0

    

        # Beams
        for beam_section_name, beam_section_tag in self.beam_section_tags.items():
            beam_data = WF_Database(beam_section_name)
            beam = I_shape(beam_data.d, beam_data.tw, beam_data.bf, beam_data.tf,   
                        Fy=Steel.Fy, E=Steel.E,
                        A=beam_data.A, 
                        Ix=beam_data.Ix,Zx=beam_data.Zx,Sx=beam_data.Sx,rx=beam_data.rx,
                        Iy=beam_data.Iy,Zy=beam_data.Zy,Sy=beam_data.Sy,ry=beam_data.ry,
                        J=beam_data.J,Cw=beam_data.Cw,rts=beam_data.rts,ho=beam_data.ho)
            beam.build_ops_fiber_section(beam_section_tag,
                                        start_material_id=Steel.mat_tag,
                                        mat_type=mat_type,
                                        nfy=self.nfy, nfx=self.nfx,
                                        frc=frc,num_regions=self.num_regions,
                                        stiffness_reduction=self.stiffness_reduction,strength_reduction=self.strength_reduction,
                                        axis='x')
            Steel.mat_tag += 2 * self.num_regions + 2
            ops.beamIntegration("Lobatto", beam_section_tag, beam_section_tag, self.nip)
            setattr(self, beam_section_name, beam)

        # Columns
        for column_section_name, (column_section_tag, axis) in self.column_section_tags.items():
            column_data = WF_Database(column_section_name)
            column = I_shape(column_data.d, column_data.tw, column_data.bf, column_data.tf,
                            Fy=Steel.Fy, E=Steel.E,
                            A=column_data.A,
                            Ix=column_data.Ix, Zx=column_data.Zx, Sx=column_data.Sx, rx=column_data.rx,
                            Iy=column_data.Iy, Zy=column_data.Zy, Sy=column_data.Sy, ry=column_data.ry,
                            J=column_data.J, Cw=column_data.Cw, rts=column_data.rts, ho=column_data.ho)
            column.build_ops_fiber_section(column_section_tag,
                                        start_material_id=Steel.mat_tag,
                                        mat_type=mat_type,
                                        nfy=self.nfy, nfx=self.nfx,
                                        frc=frc,num_regions=self.num_regions,
                                        stiffness_reduction=self.stiffness_reduction,strength_reduction=self.strength_reduction,
                                        axis=axis)
            Steel.mat_tag += 2 * self.num_regions + 2
            ops.beamIntegration("Lobatto", column_section_tag, column_section_tag, self.nip)
            setattr(self, column_section_name, column)


        ########### check whether Second order effects need to be included or not ############
        if self.Second_order_effects:    
            ops.geomTransf("PDelta", col_and_beam_TransTag)

        else:
            ops.geomTransf("Linear", col_and_beam_TransTag)
 





        # ----- Columns -----
        for i, column_ij_node in enumerate(self.column_connectivity):
            
            axis, storey = Moment_Frame_2D.find_axis_and_storey_for_columns(column_ij_node[1:])

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
            bay, storey = Moment_Frame_2D.find_bay_and_storey_for_beams(beam_ij_node[1:])

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

        self.roof_beams=[]
        for beams in self.beam_connectivity:
            if Moment_Frame_2D.nth_digit(beams[1],2,2)==self.no_of_stories:
                self.roof_beams.append(beams)

###########################################################################################################################
###########################################################################################################################
    # def return_fiber_strain_in_an_element(self,ele_tag):
    #     for i in range(self.nip):
    #         axial_strain, curvatureX, curvatureY = 0, 0, 0
    #         resp = ops.eleResponse(ele_tag,
    #                                'section', i+1,    ##integration point
    #                                'deformation') 
    #         print('line 605',resp)
    #         return resp[1]
    
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

    def return_max_of_P_M_M_interaction(self):
        P_M_M_interaction_all_elements=[]
        for ele in self.all_element_connectivity_section_and_bending_axes_detail:
            ele_tag=ele[0]
            ele_node_i=ele[1]
            ele_node_j=ele[2]
            ele_section_name=ele[3]        #### str  'W8X15'
            ele_bending_axis=ele[4]        #### str  'x'
            ele_type=ele[5]
            coords_i = ops.nodeCoord(ele_node_i)  # [x_i, y_i]
            coords_j = ops.nodeCoord(ele_node_j)  # [x_j, y_j]
            L = math.sqrt((coords_j[0] - coords_i[0])**2 + (coords_j[1] - coords_i[1])**2)

            section_obj = getattr(self, ele_section_name)   # retrieves the I_shape instance
            
            # Build a unique name for the member instance
            member_name = f"member_{ele_tag}_{ele_section_name}_{ele_bending_axis}"
            
            # Create the WideFlangeMember_AISC2023 instance
            member_obj = WideFlangeMember_AISC2022( section_obj,
                                                    Fy=section_obj.Fy,
                                                    E=section_obj.E,
                                                    L=L  
                                                )
            
            setattr(self, member_name, member_obj) 
            ### Required Strength and Available Strengths
            forces = ops.eleResponse(ele_tag, 'localForce')
            Pr=abs(forces[0])  ## required compressive strength
            Mr_i=forces[2]
            Mr_j=forces[5]
            max_Mr=max(abs(Mr_i),abs(Mr_j))
            Mrx=0
            Mry=0
            Mcx=member_obj.Mnx(Lb=0,Cb=1)  ##available flexural strength for bending about major axis
            Mcy=member_obj.Mny()   ##available flexural strength for bending about minor axis

            if ele_bending_axis=='x':
                Mrx=max_Mr
                if ele_type=='col':
                    Leff=L*self.no_of_elements_column
                if ele_type=='beam': 
                    Leff=L*self.no_of_elements_beam
                Pcc=member_obj.Pnc(Lcx=Leff,Lcy=0) ##available compressive strength 

            if ele_bending_axis=='y':
                Mry=max_Mr
                if ele_type=='col':
                    Leff=L*self.no_of_elements_column
                if ele_type=='beam': 
                    Leff=L*self.no_of_elements_beam
                Pcc=member_obj.Pnc(Lcx=0,Lcy=Leff) ##available compressive strength 

            ### Interaction Eqn H1-1a,H1-1b
            if Pr/Pcc>=0.2:
                P_M_M_interaction=Pr/Pcc+(8/9)*((Mrx/Mcx)+(Mry/Mcy))
            else:
                P_M_M_interaction=Pr/(2*Pcc)+((Mrx/Mcx)+(Mry/Mcy))
            P_M_M_interaction_all_elements.append(P_M_M_interaction)

            # print(ele_tag)
            # print(ele_section_name)
            # print(L)
            # print('Mnx',Mcx)
            # print('Mny',Mcy)
            # print('Pcc',Pcc)
            # print('Mrx',Mrx)
            # print('Mry',Mry)
            # print('Pr',Pr)
            # print(P_M_M_interaction)
            # print(P_M_M_interaction_all_elements)
        max_P_M_M_interaction_among_all_elements=max(P_M_M_interaction_all_elements)
        # print(max_P_M_M_interaction_among_all_elements)
        # print('\n')
        # print('\n')
        # print('\n')
        # print('\n')
        # input('Hello')
        return max_P_M_M_interaction_among_all_elements

    

    def add_vertical_dead_live_wall_loads(self,vertical_load_scale=1):   
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


    def add_lateral_wind_notional_loads(self, lateral_load_scale=1.0):
        ##Lateral Loads
        # Wind  Load 
        if self.wind_load_dirn.lower() == 'right':
            for node in self.axis_i_floor_nodes(1):
                ops.load(node, lateral_load_scale * self.Wind_load_floor * self.W_multiplier, 0, 0.0)
            for node in self.axis_i_roof_nodes(1):
                ops.load(node, lateral_load_scale * self.Wind_load_roof * self.W_multiplier, 0, 0.0)

        elif self.wind_load_dirn.lower() == 'left':
            for node in self.axis_i_floor_nodes(self.no_of_bays + 1):
                ops.load(node, -lateral_load_scale * self.Wind_load_floor * self.W_multiplier, 0, 0.0)
            for node in self.axis_i_roof_nodes(self.no_of_bays + 1):
                ops.load(node, -lateral_load_scale * self.Wind_load_roof * self.W_multiplier, 0, 0.0)

        # Notional Load
        if self.Notional_load:
            for node in self.axis_i_floor_nodes(1):
                ops.load(node, lateral_load_scale * self.floor_notional_load(), 0, 0.0)
                # print('line852',node,self.floor_notional_load())
            for node in self.axis_i_roof_nodes(1):
                ops.load(node, lateral_load_scale * self.roof_notional_load(), 0, 0.0)
                # print('line856',node,lateral_load_scale * self.roof_notional_load())            
                
        # opsv.plot_model()
        # opsv.plot_load()


    def run_load_controlled_analysis(self,steps = 10,plot_defo=False,display_reactions=False):

            
        """
        Runs gravity analysis.
        Note that the model should be built before
        calling this function.
        
        Keyword arguments:
        steps -- total number of analysis steps

        """

        ops.initialize()
        # Records the response of a number of nodes at every converged step
        self.main_node_tags=[tags[0] for tags in sorted(self.Main_Nodes,key=lambda x:x[2])]
        # Ensure output folder exists
        os.makedirs(self.Frame_id, exist_ok=True)
        filename=self.Frame_id+'/Gravity_Displacements.out'
        ops.recorder('Node', '-file', filename,
                    '-time','-node', *self.main_node_tags, '-dof',*[1] , 'disp')

        # plain constraint handler enforces homogeneous single point constraints
        ops.constraints('Plain')

        # RCM numberer uses the reverse Cuthill-McKee scheme to order the matrix equations
        ops.numberer('RCM')

        # Constructs a profileSPDSOE (Symmetric Positive Definite) system of equation object
        ops.system('ProfileSPD')

        # Uses the norm of the left hand side solution vector of the matrix equation to
        # determine if convergence has been reached
        ops.test('NormDispIncr', 1.0e-6, 100, 0, 2)

        # Uses the Newton-Raphson algorithm to solve the nonlinear residual equation
        ops.algorithm('Newton')

        # Uses LoadControl integrator object
        ops.integrator('LoadControl', 1/steps)

        # Constructs the Static Analysis object
        ops.analysis('Static')

        # Records the current state of the model
        ops.record()
        # Performs the analysis
        ops.analyze(steps)    
        Moment_Frame_2D.generate_clean_csv_file_from_messy_out_files(filename)

        # === Track horizontal displacement of top node in axis 1 ===
        '''
        By tracking the displacement of the top node, the displacement controlled analysis decides whether to 
        perform the analysis for left or right lateral displacement'''
        top_node_tag = self.axis_i_roof_nodes(1)[0]  # This returns the tag of top node in the frist axis
        disp_x = ops.nodeDisp(top_node_tag, 1)       # DOF 1 = horizontal
        print(f"Horizontal displacement of control node {top_node_tag} = {disp_x:.5f} ")
              
        if display_reactions:
            print("Reactions at base nodes:")
            total_rx = 0.0
            total_ry = 0.0
            total_Mz = 0.0

            for base_node in [n[0] for n in self.NODES_TO_FIX]:  # extract node tags from base node lists
                ops.reactions()
                rxn=ops.nodeReaction(base_node)
                rx, ry, rz = rxn[0], rxn[1], rxn[2]
                total_rx += rx
                total_ry += ry
                total_Mz += rz
                print(f"Node {base_node}: Rx = {rx:.4f}, Ry = {ry:.4f}, Mz = {rz:.4f}")

            print(f"\nTotal Reactions: Rx = {total_rx:.4f}, Ry = {total_ry:.4f}, Mz = {total_Mz:.4f}")


        # forces=ops.eleResponse(9,'localForce')
        # print("Forces",forces)
        print("Gravity analysis Done!")
        if plot_defo==True:
            opsv.plot_defo()
        else:
            pass
        
        return disp_x

    def run_displacement_controlled_analysis(self, target_disp=10, steps=20000, plot_defo=False,**kwargs):
        """
        Runs displacement-controlled analysis and plots load ratio (λ) vs. displacement and vertical reaction.

        Parameters:
            target_disp (float): Target horizontal displacement at control node 
            steps (int): Number of steps to reach target
            plot_defo (bool): Whether to plot deformed shape at end
        """
        incr_LCA= kwargs.get('incr_LCA', 0.02)          ######### LCA refers to Load Controlled Analysis
        num_steps_LCA= kwargs.get('num_steps_LCA', 20)            ######### LCA refers to Load Controlled Analysis
        steel_strain_limit = kwargs.get('steel_strain_limit', 0.05)
        eigenvalue_limit = kwargs.get('eigenvalue_limit', 0)
        P_M_M_interaction_limit=kwargs.get('P_M_M_interaction_limit',1)
        try_smaller_steps = kwargs.get('try_smaller_steps', True)
        control_dir=kwargs.get('control_dir','L')  # L for lateral and V for Vertical
        ops_analysis=kwargs.get('analysis','proportional_limit_point')
        lateral_load_scale=kwargs.get('lateral_load_scale',1)
        vertical_load_scale=kwargs.get('vertical_load_scale',1)


        # Initialize analysis results
        results = AnalysisResults()
        attributes = ['load_ratio','vertical_reaction','base_shear','control_node_displacement', 'control_node_displacement_absolute',
                      'lowest_eigenvalue','absolute_maximum_strain','max_P_M_M_interaction']
        
        for attr in attributes:
            setattr(results, attr, [])


        # Define function to find limit point
        def find_limit_point():
            if Moment_Frame_2D.print_ops_status:
                print(results.exit_message)
            if 'Moving to Displacement Controlled Analysis' in results.exit_message:
                return
            if  'Analysis Failed In Load Controlled Loading before entering Displacement controlled Loading' in results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            if 'Analysis Failed' in results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif 'Eigenvalue Limit' in results.exit_message:
                ind, x = find_limit_point_in_list(results.lowest_eigenvalue, eigenvalue_limit)
            elif 'Extreme Steel Fiber Strain Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.absolute_maximum_strain, steel_strain_limit)
            elif 'P_M_M interaction Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.max_P_M_M_interaction, P_M_M_interaction_limit)            
            else:
                raise Exception('Unknown limit point')
            results.maximum_load_ratio_at_limit_point = interpolate_list(results.load_ratio, ind, x)
            print('Line 884, Max Load Ratio',results.maximum_load_ratio_at_limit_point)


        if ops_analysis.lower()=='proportional_limit_point':

            fail_during_LCA=True

            ops.timeSeries('Linear', self.load_timeseries_counter)
            ops.pattern('Plain',self.load_pattern_counter, self.load_timeseries_counter)
            self.add_vertical_dead_live_wall_loads(vertical_load_scale=vertical_load_scale)
            self.add_lateral_wind_notional_loads(lateral_load_scale=lateral_load_scale)
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
                results.max_P_M_M_interaction.append(self.return_max_of_P_M_M_interaction())

            # endregion

            # Control node for lateral deflection
            if control_dir=='L':
                control_direction='lateral'
                control_node = self.axis_i_roof_nodes(1)[0]   
                control_dof = 1  # horizontal displacement
            # Control node for vertical deflection of beam node
            else:
                control_direction='vertical'
                control_node = self.bay_i_internal_roof_nodes(1)[0]  
                control_dof = 2  # vertical displacement
            ops.initialize()
            # Create output folder
            os.makedirs(self.Frame_id, exist_ok=True)

            ops.constraints('Plain')
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test('NormUnbalance', 1e-3, 10)
            ops.algorithm('Newton')
            ops.integrator('LoadControl',incr_LCA)  ## incr_LCA because, we do not want to apply the entire load during load controlled analysis.
            ops.analysis('Static')
            record()
            for i in range(num_steps_LCA):
                if Moment_Frame_2D.print_ops_status:
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

                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        find_limit_point()
                        return results,fail_during_LCA
                        # break

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    # if Moment_Frame_2D.print_ops_status:
                    #     print(f'Checking Steel Tensile Strain')
                    if results.absolute_maximum_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                        find_limit_point()
                        return results,fail_during_LCA
                        # break
                # Check for maximum PMM interaction value    
                if self.Elastic_analysis:
                    if P_M_M_interaction_limit is not None:
                        # if Moment_Frame_2D.print_ops_status:
                        #     print(f'Checking PMM Interaction')
                        if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                            results.exit_message = 'P_M_M interaction Limit Reached'
                            find_limit_point()
                            return results,fail_during_LCA


            dU = target_disp / steps
            if results.control_node_displacement[-1]<0 or control_dof==2:
                dU=-dU 

            ops.integrator('DisplacementControl', control_node, control_dof, dU)

            record()

            while True:
                fail_during_LCA=False
                ok = ops.analyze(1)
                if try_smaller_steps:
                    if ok != 0:
                        if Moment_Frame_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 10}')
                        ops.integrator('DisplacementControl',control_node, control_dof, dU / 10)
                        ok = ops.analyze(1)

                    if ok != 0:
                        if Moment_Frame_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 100}')
                        ops.integrator('DisplacementControl',control_node, control_dof, dU / 100)
                        ok = ops.analyze(1)

                    if ok != 0:
                        if Moment_Frame_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 1000}')
                        ops.integrator('DisplacementControl', control_node, control_dof, dU / 1000)
                        ok = ops.analyze(1)
                        if ok == 0:
                            dU = dU / 10
                            if Moment_Frame_2D.print_ops_status:
                                print(f'Changed the step size to: {dU}')

                    if ok != 0:
                        if Moment_Frame_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 10000}')
                        ops.integrator('DisplacementControl', control_node, control_dof, dU / 10000)
                        ok = ops.analyze(1)
                        if ok == 0:
                            dU = dU / 10
                            if Moment_Frame_2D.print_ops_status:
                                print(f'Changed the step size to: {dU / 10}')

                if ok != 0:
                    if Moment_Frame_2D.print_ops_status:
                        print('Trying ModifiedNewton')
                    ops.algorithm('ModifiedNewton')
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Moment_Frame_2D.print_ops_status:
                            print('ModifiedNewton worked')

                if ok != 0:
                    if Moment_Frame_2D.print_ops_status:
                        print('Trying KrylovNewton')
                    ops.algorithm('KrylovNewton')
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Moment_Frame_2D.print_ops_status:
                            print('KrylovNewton worked')

                if ok != 0:
                    if Moment_Frame_2D.print_ops_status:
                        print('Trying KrylovNewton and Greater Tolerance')
                    ops.algorithm('KrylovNewton')
                    ops.test('NormUnbalance', 1e-4, 10)
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Moment_Frame_2D.print_ops_status:
                            print('KrylovNewton worked')

                if ok == 0:
                    # Reset analysis options
                    ops.algorithm('Newton')
                    ops.test('NormUnbalance', 1e-3, 10)
                    ops.integrator('DisplacementControl', control_node, control_dof, dU)
                else:
                    print('Analysis Failed')
                    results.exit_message = 'Analysis Failed'
                    break


                record()

                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        break

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    # if Moment_Frame_2D.print_ops_status:
                    #     print(f'Checking Steel Tensile Strain')
                    if results.absolute_maximum_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                        break
                # Check for maximum PMM interaction value    
                if self.Elastic_analysis:
                    if P_M_M_interaction_limit is not None:
                        # if Moment_Frame_2D.print_ops_status:
                        #     print(f'Checking PMM Interaction')
                        if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                            results.exit_message = 'P_M_M interaction Limit Reached'
                            break

            find_limit_point()


        elif ops_analysis.lower()=='non_proportional_limit_point':
            
            fail_during_LCA=True

            ops.timeSeries('Linear', self.load_timeseries_counter)
            ops.pattern('Plain',self.load_pattern_counter, self.load_timeseries_counter)
            self.add_vertical_dead_live_wall_loads(vertical_load_scale=vertical_load_scale)
            # self.add_lateral_wind_notional_loads(lateral_load_scale=lateral_load_scale)
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
                results.max_P_M_M_interaction.append(self.return_max_of_P_M_M_interaction())

            # endregion

            # Control node for lateral deflection
            if control_dir=='L':
                control_direction='lateral'
                control_node = self.axis_i_roof_nodes(1)[0]   
                control_dof = 1  # horizontal displacement
            # Control node for vertical deflection of beam node
            else:
                control_direction='vertical'
                control_node = self.bay_i_internal_roof_nodes(1)[0]  
                control_dof = 2  # vertical displacement
            ops.initialize()
            # Create output folder
            os.makedirs(self.Frame_id, exist_ok=True)

            ops.constraints('Plain')
            ops.numberer('RCM')
            ops.system('UmfPack')
            ops.test('NormUnbalance', 1e-3, 10)
            ops.algorithm('Newton')
            ops.integrator('LoadControl',1/num_steps_LCA)  ## 1/num_steps_LCA because we want the entire vertical load to be applied before starting displacement controlled analysis.
            ops.analysis('Static')
            record()
            for i in range(num_steps_LCA):
                if Moment_Frame_2D.print_ops_status:
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

                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        find_limit_point()
                        return results,fail_during_LCA
                        # break

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    # if Moment_Frame_2D.print_ops_status:
                    #     print(f'Checking Steel Tensile Strain')
                    if results.absolute_maximum_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                        find_limit_point()
                        return results,fail_during_LCA
                        # break
                # Check for maximum PMM interaction value    
                if self.Elastic_analysis:
                    if P_M_M_interaction_limit is not None:
                        # if Moment_Frame_2D.print_ops_status:
                        #     print(f'Checking PMM Interaction')
                        if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                            results.exit_message = 'P_M_M interaction Limit Reached'
                            find_limit_point()
                            return results,fail_during_LCA
                            # break
            # find_limit_point()


            dU = target_disp / steps
            if results.control_node_displacement[-1]<0 or control_dof==2:
                dU=-dU 
            ops.loadConst('-time', 0.0)
            ops.timeSeries('Linear', self.load_timeseries_counter+1)
            ops.pattern('Plain',self.load_pattern_counter+1, self.load_timeseries_counter+1)
            self.add_lateral_wind_notional_loads()
            # self.add_vertical_dead_live_wall_loads()
            ops.integrator('DisplacementControl', control_node, control_dof, dU)

            record()
            i=1
            while True:
                print(f'Running Displacement Controlled Analysis {i}')
                print(lateral_load_scale)
                i=i+1
                fail_during_LCA=False
                ok = ops.analyze(1)
                if try_smaller_steps:
                    if ok != 0:
                        if Moment_Frame_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 10}')
                        ops.integrator('DisplacementControl',control_node, control_dof, dU / 10)
                        ok = ops.analyze(1)

                    if ok != 0:
                        if Moment_Frame_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 100}')
                        ops.integrator('DisplacementControl',control_node, control_dof, dU / 100)
                        ok = ops.analyze(1)

                    if ok != 0:
                        if Moment_Frame_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 1000}')
                        ops.integrator('DisplacementControl', control_node, control_dof, dU / 1000)
                        ok = ops.analyze(1)
                        if ok == 0:
                            dU = dU / 10
                            if Moment_Frame_2D.print_ops_status:
                                print(f'Changed the step size to: {dU}')

                    if ok != 0:
                        if Moment_Frame_2D.print_ops_status:
                            print(f'Trying the step size of: {dU / 10000}')
                        ops.integrator('DisplacementControl', control_node, control_dof, dU / 10000)
                        ok = ops.analyze(1)
                        if ok == 0:
                            dU = dU / 10
                            if Moment_Frame_2D.print_ops_status:
                                print(f'Changed the step size to: {dU / 10}')

                if ok != 0:
                    if Moment_Frame_2D.print_ops_status:
                        print('Trying ModifiedNewton')
                    ops.algorithm('ModifiedNewton')
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Moment_Frame_2D.print_ops_status:
                            print('ModifiedNewton worked')

                if ok != 0:
                    if Moment_Frame_2D.print_ops_status:
                        print('Trying KrylovNewton')
                    ops.algorithm('KrylovNewton')
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Moment_Frame_2D.print_ops_status:
                            print('KrylovNewton worked')

                if ok != 0:
                    if Moment_Frame_2D.print_ops_status:
                        print('Trying KrylovNewton and Greater Tolerance')
                    ops.algorithm('KrylovNewton')
                    ops.test('NormUnbalance', 1e-4, 10)
                    ok = ops.analyze(1)
                    if ok == 0:
                        if Moment_Frame_2D.print_ops_status:
                            print('KrylovNewton worked')

                if ok == 0:
                    # Reset analysis options
                    ops.algorithm('Newton')
                    ops.test('NormUnbalance', 1e-3, 10)
                    ops.integrator('DisplacementControl', control_node, control_dof, dU)
                else:
                    print('Analysis Failed')
                    results.exit_message = 'Analysis Failed'
                    break


                record()

                # Check for lowest eigenvalue less than zero
                if eigenvalue_limit is not None:
                    if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                        results.exit_message = 'Eigenvalue Limit Reached'
                        break

                # Check for strain in extreme steel fiber
                if steel_strain_limit is not None:
                    # if Moment_Frame_2D.print_ops_status:
                    #     print(f'Checking Steel Tensile Strain')
                    if results.absolute_maximum_strain[-1] > steel_strain_limit:
                        results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                        break
                # Check for maximum PMM interaction value    
                if self.Elastic_analysis:
                    if P_M_M_interaction_limit is not None:
                        # if Moment_Frame_2D.print_ops_status:
                        #     print(f'Checking PMM Interaction')
                        if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                            results.exit_message = 'P_M_M interaction Limit Reached'
                            break

            find_limit_point()
       
       
        else:
            raise Exception('Give valid ops_analysis option')
        # Optional: plot deformed shape
        if plot_defo:
            try:
                import opsvis
                opsvis.plot_defo()
            except:
                print("opsvis not available for deformation plotting.")
        results.control_node_displacement_absolute[:] = [abs(x) if x is not None else None
                                        for x in results.control_node_displacement]




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
        plot_undeformed_2d()

    def display_node_coords(self):
        get_node_coords_and_disp()


    