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


    def __init__(self,width_of_bay,storey_height, no_of_elements_column, no_of_elements_beam):

        self.bay_width=[0]+width_of_bay
        self.storey_height = [0] + storey_height
        print(self.bay_width)
        # Assigning values to instance variables
        self.no_of_bays = len(self.bay_width)-1
        self.no_of_stories = len(self.storey_height)-1
        self.no_of_elements_column = no_of_elements_column
        self.no_of_elements_beam = no_of_elements_beam
        self.no_of_nodes_column=self.no_of_elements_column-1
        self.no_of_nodes_beam=self.no_of_elements_beam-1
        self.node_to_node_height_column=[storey_height/self.no_of_elements_column for storey_height in self.storey_height]
        self.node_to_node_height_column=self.node_to_node_height_column[1:]
        self.node_to_node_length_beam=[bay_width/self.no_of_elements_beam for bay_width in self.bay_width]
        self.node_to_node_length_beam=self.node_to_node_length_beam[1:]
        print(self.node_to_node_height_column)
        print(self.node_to_node_length_beam)
        self.Main_Nodes=[]

        # Input validation
        if not (0<= self.no_of_bays < 10):
            raise ValueError("Number of bays must be from 1 to 9.")

        if not(1<=self.no_of_stories < 100):
            raise ValueError("Number of stories must be from 1 to 99.")

        if not(1<=self.no_of_elements_column < 10):
            raise ValueError("Number of elements in columns must be from 1 to 9.")

        if not(1<=self.no_of_elements_beam < 9):
            raise ValueError("Number of elements in beams must be from 1 to 9.")


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

    def build_ops_model(self):
        nip = 7
        ops.wipe()
        
        ops.model('basic','-ndm',2,'-ndf',3)
        
        for single_node in self.all_nodes:
            ops.node(single_node[0],single_node[1],single_node[2])
    

        for base_nodes in self.NODES_TO_FIX:
            ops.fix(base_nodes[0],1,1,1)

        Reinf_steel={'Rebar_steel_tag':3,'fy':415*Mpa,'Es':200000*Mpa,'b':0.002,'R0':15,'cR1':0.925,'cR2':0.15}
        ops.uniaxialMaterial("Steel02",Reinf_steel['Rebar_steel_tag'],Reinf_steel['fy'],Reinf_steel['Es'],
                            Reinf_steel['b'],Reinf_steel['R0'],Reinf_steel['cR1'],Reinf_steel['cR2'])
        colSecTag= 1
        # # d=13.8*inch
        # d=section.wide_flange_database['W14X48']['d']*inch
        # # tw=0.34*inch
        # tw=section.wide_flange_database['W14X48']['tw']*inch
        # bf=section.wide_flange_database['W14X48']['bf']*inch
        # tf=section.wide_flange_database['W14X48']['tf']*inch
        # # A=14.1*inch*inch
        # A=section.wide_flange_database['W14X48']['A']*(inch**2)
        # # I=484*inch*inch*inch*inch
        # I1=section.wide_flange_database['W14X48']['Ix']*(inch**4)
        # I2=section.wide_flange_database['W14X48']['Iy']*(inch**4)
        # E=200000*Mpa

        # G=77221*Mpa
        # alphaY=(d*tw)/A

        # ops.section('Elastic',colSecTag,E,A,I1,G,alphaY)
        W14x48=WF_section('W14X48')
        #  sec_tag, mat_tag, nf_dw, nf_tw, nf_bf, nf_tf,section_name,plot=True
        W14x48.W_Fiber_section_minor_axis(colSecTag,Reinf_steel['Rebar_steel_tag'], *[4, 2, 4, 2],'Steel_Section',plot=True)

        colTransTag = 1
        ops.geomTransf("PDelta", colTransTag)
        beam_integration_tag=1
        ops.beamIntegration("Lobatto",beam_integration_tag , colSecTag, nip)

        for column_ij_node in self.column_connectivity:
            ops.element('forceBeamColumn',column_ij_node[0],*column_ij_node[1:],colTransTag,beam_integration_tag)

        for beam_ij_node in self.beam_connectivity:
            ops.element('forceBeamColumn',beam_ij_node[0],*beam_ij_node[1:],colTransTag,beam_integration_tag)

    def plot_model(self):
        plot_undeformed_2d()

    def display_node_coords(self):
        get_node_coords_and_disp()



  
#width_of_bay,storey_height, no_of_elements_column, no_of_elements_beam
# Frame=Moment_Frame_2D(3,3,1,3.5,3,4)
Frame=Moment_Frame_2D([3,2,3,6],[3,3.5,3,5],3,4)
Frame.generate_Nodes_and_Element_Connectivity()
# print(Frame.Main_Nodes)
# print(Frame.NODES_TO_FIX)
# print(Frame.column_intermediate_nodes)
# print(Frame.column_connectivity)
# print(Frame.beam_intermediate_nodes)
# print(Frame.beam_connectivity)
# print(Frame.all_nodes)

Frame.build_ops_model()
Frame.plot_model()
Frame.display_node_coords()
    