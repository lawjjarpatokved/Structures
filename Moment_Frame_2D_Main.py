import numpy as np
import math as math
import pandas as pd
import openseespy.opensees as ops
import time
import matplotlib.pyplot as plt
import csv
import os
import scipy as scp
import opsvis as opsv
import PyPonding.structures.wide_flange as wf
import libdenavit.section.database.aisc as section
from Units import plot_undeformed_2d,get_node_coords_and_disp

######## UNITS ##################################
m=1.0; mm=0.001*m; cm=0.01*m; km=1000.0*m; inch= 0.0254*m;ft=0.3048*m
KN=1.0; N=0.001*KN
Mpa= 10**3* KN/m**2; Gpa= 10**6* KN/m**2
sec=1
tonne=1*KN*(sec**2)/m; kg=0.001*tonne
#################################################
density_of_steel=7850*kg/(m**3)
g=9.81*m/(sec**2)
E=200000*Mpa
G=77221*Mpa


class WF_section:

    def __init__(self,Section_name):

        self.section=Section_name
        self.d=section.wide_flange_database[self.section]['d']*inch
        self.tw=section.wide_flange_database[self.section]['tw']*inch
        self.bf=section.wide_flange_database[self.section]['bf']*inch
        self.tf=section.wide_flange_database[self.section]['tf']*inch
        self.A=section.wide_flange_database[self.section]['A']*(inch**2)
        self.I1=section.wide_flange_database[self.section]['Ix']*(inch**4)
        self.I2=section.wide_flange_database[self.section]['Iy']*(inch**4)

    def W_Fiber_section_minor_axis(self, sec_tag, mat_tag, nf_dw, nf_tw, nf_bf, nf_tf,section_name,plot=True):
        
        """
        Creates W-Section based on nominal dimension and generates
        fibers over it
        TCL version by: Remo M. de Souza, 08/99
        Python version by: Amir Hossein Namadchi, 2022
        
        Keyword arguments:
        section -- a dict type containing section info
        sec_tag -- Section tag
        mat_tag -- Material Tag
        nf_dw -- Number of fibers along web depth 
        nf_tw -- Number of fibers along web thickness
        nf_bf -- Number of fibers along flange width
        nf_tf -- Number of fibers along flange thickness    
        """
        
        dw = self.d - 2*self.tf
        y1, y2, y3, y4 = -self.d/2, -dw/2, dw/2,self.d/2
        z1, z2, z3, z4 = -self.bf/2, -self.tw/2, self.tw/2, self.bf/2
        
        ops.section('Fiber', sec_tag)
        ops.patch('quad', mat_tag, nf_bf, nf_tf,
            *[y1,z4], *[y1,z1], *[y2,z1], *[y2,z4])
        ops.patch('quad', mat_tag, nf_tw, nf_dw,
            *[y2,z3], *[y2,z2], *[y3,z2], *[y3,z3])    
        ops.patch('quad', mat_tag, nf_bf, nf_tf,
            *[y3,z4], *[y3,z1], *[y4,z1], *[y4,z4])
        
        if plot:
            fib_sec = [

                    ['section', 'Fiber', sec_tag, '-GJ', 1.0e4],

                    ["patch","quad", mat_tag, nf_bf, nf_tf,*[y1,z4], *[y1,z1], *[y2,z1], *[y2,z4]],
                    ["patch","quad",mat_tag, nf_tw, nf_dw,*[y2,z3], *[y2,z2], *[y3,z2], *[y3,z3]],
                    ["patch","quad", mat_tag, nf_bf, nf_tf,*[y3,z4], *[y3,z1], *[y4,z1], *[y4,z4]]

                    ]
            
            matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
            opsv.plot_fiber_section(fib_sec, matcolor=matcolor)
            plt.axis('equal')
            plt.title(section_name)
            plt.show()


class Moment_Frame_2D:


    def __init__(self,width_of_bay,storey_height,
                no_of_elements_column, no_of_elements_beam,
                 beam_section,column_section,load_combination_multipliers,
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
        self.make_beam_section_detail_uniform()
        self.make_column_section_detail_uniform()

        self.kwargs=kwargs

        defaults={'support':'All_Fixed',
                  'nip':4,
                  'D_floor_intensity':0,
                  'D_roof_intensity':0,
                  'L_floor_intensity':0,
                  'L_roof_intensity':0,
                  'Inelastic_analysis':True,
                  'Second_order_effects':True}
        
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
    def nth_digit(num, n):
        """
        Returns the nth digit (from left, 1-based index) of a number.
        For example, nth_digit(1234, 2) returns 2.
        """
        num_str = str(abs(num))  # Make sure it's positive and convert to string
        if n <= 0 or n > len(num_str):
            raise ValueError(f"n = {n} is out of bounds for number {num}")
        return int(num_str[n - 1])

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
        storey = Moment_Frame_2D.nth_digit(first_node, 3)
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
        storey = Moment_Frame_2D.nth_digit(first_node, 3) + 1  
        ''' 1 is added because my assumption is that column in storey
          1 means the column that support the 1st storey slab
        '''
        return axis, storey


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
       


        self.all_nodes=self.Main_Nodes+self.column_intermediate_nodes+self.beam_intermediate_nodes


    

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

        # Assign unique tags to unique section names
        unique_sections = sorted(set(nested_dict.values()))
        self.column_section_tags = {section: idx + 1+self.no_of_beam_sections for idx, section in enumerate(unique_sections)}

        # Update the nested dictionary to include the tag
        updated_nested_dict = {
            key: (value, self.column_section_tags[value]) for key, value in nested_dict.items()
        }

        # Update the main beam_section
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


        Reinf_steel={'Rebar_steel_tag':3,'fy':415*Mpa,'Es':200000*Mpa,'b':0.002,'R0':15,'cR1':0.925,'cR2':0.15}
        ops.uniaxialMaterial("Steel02",Reinf_steel['Rebar_steel_tag'],Reinf_steel['fy'],Reinf_steel['Es'],
                            Reinf_steel['b'],Reinf_steel['R0'],Reinf_steel['cR1'],Reinf_steel['cR2'])
        


        col_and_beam_TransTag = 1
        beam_integration_tag=1
        column_integration_tag=2


        if not self.Inelastic_analysis:
            for beam_section_name,beam_section_tag in self.beam_section_tags.items():
                beam=WF_section(beam_section_name)
                A=beam.A
                I=beam.I1 ### This is for major axis bending, for minor axis use beam.I2
                d=beam.d
                tw=beam.tw
                alphaY=(d*tw)/A
                ops.section('Elastic',beam_section_tag,E,A,I,G,alphaY)
                ops.beamIntegration("Lobatto",beam_section_tag , beam_section_tag, self.nip)

            for column_section_name, column_section_tag in self.column_section_tags.items():
                column = WF_section(column_section_name)
                A = column.A
                I = column.I1  ### This is for major axis bending, for minor axis use column.I2
                d = column.d
                tw = column.tw
                alphaY = (d * tw) / A
                ops.section('Elastic', column_section_tag, E, A, I, G, alphaY)
                ops.beamIntegration("Lobatto",column_section_tag , column_section_tag, self.nip)


        elif self.Inelastic_analysis:

            for beam_section_name,beam_section_tag in self.beam_section_tags.items():
                beam=WF_section(beam_section_name)
                                                
                beam.W_Fiber_section_minor_axis(beam_section_tag,               #  sec_tag
                                                Reinf_steel['Rebar_steel_tag'], # mat_tag
                                                *[4, 2, 4, 2],                  #  nf_dw, nf_tw, nf_bf, nf_tf
                                                'Steel_Section',                # section_name
                                                plot=True)                      # plot=True
                
                ops.beamIntegration("Lobatto",beam_section_tag , beam_section_tag, self.nip)


            for column_section_name,column_section_tag in self.column_section_tags.items():
                column=WF_section(column_section_name)
                                                
                column.W_Fiber_section_minor_axis(column_section_tag,           #  sec_tag
                                                Reinf_steel['Rebar_steel_tag'], # mat_tag
                                                *[4, 2, 4, 2],                  #  nf_dw, nf_tw, nf_bf, nf_tf
                                                'Steel_Section',                # section_name
                                                plot=True)                      # plot=True
                
                ops.beamIntegration("Lobatto",column_section_tag , column_section_tag, self.nip)

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

            # Append section name to beam entry
            self.beam_connectivity[i] = beam_ij_node + [section_name]

        self.roof_beams=[]
        for beams in self.beam_connectivity:
            if Moment_Frame_2D.nth_digit(beams[1],3)==self.no_of_stories:
                self.roof_beams.append(beams)

        #### This part of code adds the self weight of all beams and columns, also any additional dead or live load 
        load_timeseries_counter=1
        load_pattern_counter=1
        ops.timeSeries('Linear',load_timeseries_counter)
        ops.pattern('Plain',load_pattern_counter,load_timeseries_counter)
        load_timeseries_counter+=1
        load_pattern_counter+=1

        for columns in self.column_connectivity:
            fiber_section=WF_section(columns[3])
            Area=fiber_section.A
            ops.eleLoad('-ele',columns[0],'-type','-beamUniform',0,-Area*density_of_steel*g)


        for beams in self.beam_connectivity:
            fiber_section=WF_section(beams[3])
            Area=fiber_section.A
            # ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-Area*density_of_steel*g,0)
            if beams not in self.roof_beams:
                ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-self.D_multiplier*self.D_floor_intensity,0) 
                ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-self.L_multiplier*self.L_floor_intensity,0)
            else:
                ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-self.D_multiplier*self.D_roof_intensity,0)
                ops.eleLoad('-ele',beams[0],'-type','-beamUniform',-self.L_r_multiplier*self.L_roof_intensity,0)

            

        


        opsv.plot_model()
        opsv.plot_load(nep=2)

    # def apply_dead_and_live_loads(self):



    def plot_model(self):
        plot_undeformed_2d()

    def display_node_coords(self):
        get_node_coords_and_disp()


    