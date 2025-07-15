from Moment_Frame_2D_Main import *

beam_section={'common_and_exceptions':{'common':'W14X48',
                                        '(2,1)':'W14X43',
                                        '(2,2)':'W14X53'}}

# beam_section={'same_for_storey':{'1':'W14X48',
#                                  '2':'W14X43',
#                                  '3':'W14X109',
#                                  '4':'W14X53'}}

column_section={'common_and_exceptions':{'common':'W14X159',
                                        '(2,1)':'W14X61'}}

# column_section={'same_for_storey':{'1':'W14X48',
#                                  '2':'W14X43',
#                                  '3':'W14X109',
#                                  '4':'W14X53'}}

support='fpfpp'     #### f means fixed and p means pinned. Pass it as a string.
# load_comb_multipliers=[1.2,1.6,0.5,0]            #### [D,L,Lr,W]
load_comb_multipliers=[1,1,1,1]
D_floor_intensity=2                              #### if no unit is used , it is kN/m
D_roof_intensity=3
L_floor_intensity=1.5                        
L_roof_intensity=2.5


Frame=Moment_Frame_2D([3,2,3,6],[3,3.5,3,5],2,2,
                      beam_section=beam_section,
                      column_section=column_section,nip=3,
                      D_floor_intensity=D_floor_intensity,D_roof_intensity=D_roof_intensity,
                      L_floor_intensity=L_floor_intensity,L_roof_intensity=L_roof_intensity,
                      load_combination_multipliers=load_comb_multipliers)

Frame.generate_Nodes_and_Element_Connectivity()

# print(Frame.Main_Nodes)
# print(Frame.NODES_TO_FIX)
# print(Frame.column_intermediate_nodes)
# print(Frame.column_connectivity)
# print(Frame.beam_intermediate_nodes)
# print(Frame.beam_connectivity)
# print(Frame.beam_section_tags)
# print(Frame.beam_section)
# print(Frame.beam_case)
# print(Frame.no_of_beam_sections)
# print(Frame.nip)
# print(Frame.column_section_tags)
# print(Frame.column_section)
# print(Frame.column_case)
# print(Frame.no_of_column_sections)
# print(Frame.all_nodes)

Frame.build_ops_model()
print(Frame.roof_beams)
print(Frame.column_connectivity)

# print(Frame.beam_connectivity)
Frame.plot_model()
# Frame.display_node_coords()