from Moment_Frame_2D_Main import *
from Ziemian_database import Frame
from Ziemian_database import Frame_data
Frame_number=13
dict=Frame[str(Frame_number)]
Frame=Frame_data(dict)

geometric_imperfection_ratio=1/500
if geometric_imperfection_ratio>0:
    wind_load_dirn='right'
else:
    wind_load_dirn='left'


Frame=Moment_Frame_2D(Frame.bay_width,Frame.story_height,Frame.column_no_of_ele,Frame.beam_no_of_ele,
                      beam_section=Frame.beam_section,
                      column_section=Frame.column_section,support=Frame.support,
                      D_floor_intensity=Frame.D_floor_intensity,D_roof_intensity=Frame.D_roof_intensity,
                      L_floor_intensity=Frame.L_floor_intensity,L_roof_intensity=Frame.L_roof_intensity,
                      Wind_load_floor=Frame.Wind_load_floor,Wind_load_roof=Frame.Wind_load_roof,
                      Wall_load=Frame.Wall_load,
                      load_combination_multipliers=Frame.load_comb_multipliers,Frame_id=Frame.Frame_id,
                      Residual_Stress=True,Inelastic_analysis=False,Second_order_effects=False,
                      stiffness_reduction=0.8,nip=10,
                      wind_load_dirn=wind_load_dirn
                      )

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
# print(Frame.beam_nodes)
# print(Frame.axis_i_floor_nodes(2))
# print(Frame.bay_i_internal_floor_nodes(2))
# print(Frame.bay_i_internal_roof_nodes(2))
# print(Frame.bay_i_floor_load(1))
# print(Frame.bay_i_roof_load(2))



# Frame.plot_model()
Frame.create_distorted_nodes_and_element_connectivity(geometric_imperfection_ratio=geometric_imperfection_ratio)
# print(Frame.Main_Nodes)
Frame.build_ops_model()
# Frame.plot_model()
# print(Frame.roof_beams)
# print(Frame.column_connectivity)
Frame.run_gravity_analysis(steps=1000)
# print(Frame.beam_connectivity)
Frame.plot_model()
# Frame.display_node_coords()