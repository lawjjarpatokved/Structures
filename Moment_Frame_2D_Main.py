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

#################################################
density_of_steel=7850*kg/(m**3)
g=9.81*m/(sec**2)
E=29000*ksi
G=77221*Mpa
fy=36*ksi*0.9
Hk = 0.001*E           # Kinematic hardening modulus

class Steel_Material:
    def __init__(self,mat_tag,E,fy,G,Hk,density):
        self.mat_tag=mat_tag
        self.E=E
        self.fy=fy
        self.G=G
        self.Hk=Hk
        self.density=density
        self.b=Hk / (E + Hk)

class Moment_Frame_2D:


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
        self.beam_section=beam_section
        self.column_section=column_section   
        self.Main_Nodes=[]
        self.D_multiplier=load_combination_multipliers[0]      ### Dead Load multiplier
        self.L_multiplier=load_combination_multipliers[1]      ### Live Load multiplier
        self.L_r_multiplier=load_combination_multipliers[2]    ### Roof Live Load multiplier
        self.W_multiplier=load_combination_multipliers[3]      ### Wind Load multiplier
        self.Frame_id=Frame_id
        self.make_beam_section_detail_uniform()
        self.make_column_section_detail_uniform()

        self.kwargs=kwargs

        defaults={'support':'All_Fixed',
                  'nip':4,
                  'D_floor_intensity':0,
                  'D_roof_intensity':0,
                  'L_floor_intensity':0,
                  'L_roof_intensity':0,
                  'Wind_load_floor':0,
                  'Wind_load_roof':0,
                  'Wall_load':0,
                  'Inelastic_analysis':True,
                  'Second_order_effects':True,
                  'Residual_Stress':True,
                  'stiffness_reduction':1,
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
        # Read file (whitespace-delimited)
        df = pd.read_csv(path, sep=r'\s+', header=None)

        # Select the last row, transpose, and convert to inches
        df = (df.iloc[-1] / 0.0254).to_frame()

        # Reset index and rename column for clarity
        df.reset_index(drop=True, inplace=True)
        df.columns = ['Value_in_inches']

        # Build output path in the same folder
        folder = os.path.dirname(path)
        base_name = os.path.splitext(os.path.basename(path))[0]
        output_path = os.path.join(folder, base_name + '_clean.csv')

        # Save as clean CSV
        df.to_csv(output_path, index=False)


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

    def create_distorted_nodes_and_element_connectivity(self,geometric_imperfection_ratio=None):
        ratio = geometric_imperfection_ratio if geometric_imperfection_ratio is not None else self.geometric_imperfection_ratio
        for i in range(len(self.all_nodes)):
                self.all_nodes[i][1] = self.all_nodes[i][1] +ratio * self.all_nodes[i][2]

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
    


        for base_nodes, support_condition in zip(self.NODES_TO_FIX, self.support_condition):
            node_tag = base_nodes[0]
            if support_condition.upper() == 'F':       # Fixed: UX, UY, RZ all fixed
                ops.fix(node_tag, 1, 1, 1)
            elif support_condition.upper() == 'P':     # Pinned: UX and UY fixed, RZ free
                ops.fix(node_tag, 1, 1, 0)
            else:                                      # Wrong condition
                raise ValueError(f"Unsupported support condition{support_condition}.Expected 'F' or 'P'.")


        # Reinf_steel={'Rebar_steel_tag':1,'fy':Fy,'Es':E,'b':0.002,'R0':15,'cR1':0.925,'cR2':0.15}
        # ops.uniaxialMaterial("Steel02",Reinf_steel['Rebar_steel_tag'],Reinf_steel['fy'],Reinf_steel['Es'],
        #                     Reinf_steel['b'],Reinf_steel['R0'],Reinf_steel['cR1'],Reinf_steel['cR2'])
        Steel=Steel_Material(1,E=E,fy=fy,G=G,Hk=Hk,density=density_of_steel)
        # ops.uniaxialMaterial('Steel01', Steel.mat_tag, Steel.Fy, Steel.E, Steel.b)
        


        col_and_beam_TransTag = 1



        if not self.Inelastic_analysis:
            for beam_section_name,beam_section_tag in self.beam_section_tags.items():
                beam_=WF_Database(beam_section_name)
                beam=I_shape(beam_.d,beam_.tw,beam_.bf,beam_.bf,fy=Steel.fy,E=Steel.E,Hk=Steel.Hk,A=beam_.A,Ix=beam_.Ix,Iy=beam_.Iy)
                A=beam.A
                I=beam.Ix  
                d=beam.d
                tw=beam.tw
                alphaY=(d*tw)/A
                # print(beam_section_name, A, I*2402509.61)
                ops.section('Elastic',beam_section_tag,beam.E*self.stiffness_reduction,A,I,Steel.G,alphaY)
                ops.beamIntegration("Lobatto",beam_section_tag , beam_section_tag, self.nip)

            for column_section_name, (column_section_tag, axis) in self.column_section_tags.items():
                column_ = WF_Database(column_section_name)  # Fetch geometry
                column = I_shape(column_.d, column_.tw, column_.bf, column_.tf,
                                fy=Steel.fy,E=Steel.E,Hk=Steel.Hk,
                                A=column_.A, Ix=column_.Ix, Iy=column_.Iy)  # Use the I_shape class

                A = column.A
                I = column.Ix if axis == 'x' else column.Iy
                d = column.d
                tw = column.tw
                alphaY = (d * tw) / A

                ops.section('Elastic', column_section_tag, column.E * self.stiffness_reduction, A, I, Steel.G, alphaY)
                ops.beamIntegration("Lobatto", column_section_tag, column_section_tag, self.nip)


        elif self.Inelastic_analysis:
            if self.Residual_Stress:
                frc=-0.3*Steel.fy
            else:
                frc=0

            for beam_section_name, beam_section_tag in self.beam_section_tags.items():
                beam_data = WF_Database(beam_section_name)
                beam = I_shape(beam_data.d, beam_data.tw, beam_data.bf, beam_data.tf,
                               fy=Steel.fy,E=Steel.E,Hk=Steel.Hk,
                            A=beam_data.A, Ix=beam_data.Ix, Iy=beam_data.Iy)
                
                # beam= I_shape.from_database(beam_section_name,fy=Steel.fy,E=Steel.E,Hk=Steel.Hk)
                
                beam.build_ops_fiber_section(beam_section_tag,        # sec_tag
                                                Steel.mat_tag,           # material       # section_name
                                                mat_type='Steel01',
                                                nfy=20,nfx=20,
                                                frc=frc,
                                                axis='x'
                                                )
    
                
                Steel.mat_tag += 2 * beam.num_regions + 2  # Avoid material tag overlap
                ops.beamIntegration("Lobatto", beam_section_tag, beam_section_tag, self.nip)

            for column_section_name, (column_section_tag, axis) in self.column_section_tags.items():
                column_data = WF_Database(column_section_name)
                column = I_shape(column_data.d, column_data.tw, column_data.bf, column_data.tf,
                               fy=Steel.fy,E=Steel.E,Hk=Steel.Hk,
                                A=column_data.A, Ix=column_data.Ix, Iy=column_data.Iy)

                if axis == 'x':
                    column.build_ops_fiber_section(column_section_tag,
                                                Steel.mat_tag,           # material       # section_name
                                                mat_type='Steel01',
                                                nfy=20,nfx=20,
                                                frc=frc,
                                                axis='x')
                else:
                    column.build_ops_fiber_section(column_section_tag,
                                                Steel.mat_tag,           # material       # section_name
                                                mat_type='Steel01',
                                                nfy=20,nfx=20,
                                                frc=frc,
                                                axis='y')

                Steel.mat_tag += 2 * column.num_regions + 2
                ops.beamIntegration("Lobatto", column_section_tag, column_section_tag, self.nip)
        else:
            raise ValueError('Please give correct input for Inelastic_analysis')



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
                else:
                    # print(f"node {column_ij_node[1]} and {column_ij_node[2]} are not in axis {axis} and storey {storey}, which is common.")
                    section_name = self.column_section[self.column_case]['common'][0]
                    section_tag  = self.column_section[self.column_case]['common'][1]

            elif self.column_case == 'same_for_storey':
                key = str(storey)
                if key in self.column_section[self.column_case]:
                    section_name = self.column_section[self.column_case][key][0]
                    section_tag  = self.column_section[self.column_case][key][1]
                else:
                    raise KeyError(f"No column section defined for storey {storey}")

            else:
                raise ValueError(f"Unsupported column_case: {self.column_case}. Possible issues: lowercase, spelling error.")

            # Create column element
            eleTag = column_ij_node[0]
            ops.element('forceBeamColumn', eleTag, *column_ij_node[1:], col_and_beam_TransTag, section_tag)
            # print(f"Defined element {eleTag} between nodes {column_ij_node[1]} and {column_ij_node[2]}")
            # Append section name to column entry
            self.column_connectivity[i] = column_ij_node + [section_name]





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
            ops.element('forceBeamColumn', eleTag, *beam_ij_node[1:], col_and_beam_TransTag, section_tag)
            # print(f"Defined element {eleTag} between nodes {column_ij_node[1]} and {column_ij_node[2]}")

            # Append section name to beam entry
            self.beam_connectivity[i] = beam_ij_node + [section_name]

        self.roof_beams=[]
        for beams in self.beam_connectivity:
            if Moment_Frame_2D.nth_digit(beams[1],2,2)==self.no_of_stories:
                self.roof_beams.append(beams)

###########################################################################################################################
###########################################################################################################################
        #### This part of code adds the self weight of all beams and columns, also any additional dead or live load 
        load_timeseries_counter=1
        load_pattern_counter=1
        ops.timeSeries('Linear',load_timeseries_counter)
        ops.pattern('Plain',load_pattern_counter,load_timeseries_counter)
        load_timeseries_counter+=1
        load_pattern_counter+=1

        # for columns in self.column_connectivity:
        #     fiber_section=I_shape(columns[3],Steel.E,Steel.Fy,Steel.Hk)
        #     Area=fiber_section.A
        #     ops.eleLoad('-ele',columns[0],'-type','-beamUniform',0,-Area*density_of_steel*g)


        # for beams in self.beam_connectivity:
        #     fiber_section=I_shape(beams[3],Steel.E,Steel.Fy,Steel.Hk)
        #     Area=fiber_section.A
        #     ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-Area*density_of_steel*g,0)
        #     if beams not in self.roof_beams:
        #         ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-self.D_multiplier*self.D_floor_intensity,0) 
        #         ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-self.L_multiplier*self.L_floor_intensity,0)
        #     else:
        #         ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-self.D_multiplier*self.D_roof_intensity,0)
        #         ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-self.L_r_multiplier*self.L_roof_intensity,0)


########## Dead and Live Load #########################################
        for bay in range(1, self.no_of_bays + 1):
            loaded_nodes_floor = self.bay_i_internal_floor_nodes(i=bay)
            # print(loaded_nodes_floor)
            load_value_floor = self.bay_i_floor_load(i=bay)
            # print(load_value_floor)
            for node in loaded_nodes_floor:
                ops.load(node, 0.0, load_value_floor[0], 0.0)



            loaded_nodes_roof = self.bay_i_internal_roof_nodes(bay)
            load_value_roof = self.bay_i_roof_load(bay)
            for node in loaded_nodes_roof:
                ops.load(node, 0.0, load_value_roof[0], 0.0)


        for axis in range(1, self.no_of_bays + 2):
            loaded_nodes_floor = self.axis_i_floor_nodes(axis)
            load_value_floor = self.axis_i_floor_load(axis)
            for node in loaded_nodes_floor:
                ops.load(node, 0.0, load_value_floor[0], 0.0)

            loaded_nodes_roof = self.axis_i_roof_nodes(axis)
            load_value_roof = self.axis_i_roof_load(axis)
            for node in loaded_nodes_roof:
                ops.load(node, 0.0, load_value_roof[0], 0.0)



########## Wind or Lateral Load #########################################
        if self.wind_load_dirn.lower()=='right':
            for node in self.axis_i_floor_nodes(1):
                    # print(node,self.Wind_load_floor*self.W_multiplier)
                    ops.load(node, self.Wind_load_floor*self.W_multiplier, 0, 0.0)

        elif self.wind_load_dirn.lower()=='left':
            for node in self.axis_i_floor_nodes(self.no_of_bays+1):

                    ops.load(node, -self.Wind_load_floor*self.W_multiplier, 0, 0.0)


        if self.wind_load_dirn.lower()=='right':
            for node in self.axis_i_roof_nodes(1):
                    # print(node,self.Wind_load_roof*self.W_multiplier)
                    ops.load(node, self.Wind_load_roof*self.W_multiplier, 0, 0.0)

        elif self.wind_load_dirn.lower()=='left':
            for node in self.axis_i_roof_nodes(self.no_of_bays+1):

                    ops.load(node, -self.Wind_load_roof*self.W_multiplier, 0, 0.0)

############# Wall Load ############################################
        for node in self.axis_i_floor_nodes(1):
            wall_load=self.Wall_load*self.D_multiplier
            ops.load(node,0,-wall_load,0)

        for node in self.axis_i_floor_nodes(self.no_of_bays+1):
            wall_load=self.Wall_load*self.D_multiplier
            ops.load(node,0,-wall_load,0)

        for node in self.axis_i_roof_nodes(1):
            wall_load=self.Wall_load*self.D_multiplier
            ops.load(node,0,-wall_load/2,0)

        for node in self.axis_i_roof_nodes(self.no_of_bays+1):
            wall_load=self.Wall_load*self.D_multiplier
            ops.load(node,0,-wall_load/2,0)




        opsv.plot_model()
        # opsv.plot_load()

    def plot_all_fiber_section_in_the_model(self):
        for sec_tag in self.beam_section_tags.values():
            I_shape.plot_fiber_section(section_id=sec_tag)

        for sec_tag,_ in self.column_section_tags.values():
            I_shape.plot_fiber_section(section_id=sec_tag)



    def run_gravity_analysis(self,steps = 10,plot_defo=False):
            
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


        forces=ops.eleResponse(9,'localForce')
        print("Forces",forces)
        print("Gravity analysis Done!")
        if plot_defo==True:
            opsv.plot_defo()
        else:
            pass

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


    