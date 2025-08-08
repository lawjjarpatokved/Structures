from Moment_Frame_2D_Main import *
from Ziemian_database import Frame
from Ziemian_database import Frame_data
from libdenavit.OpenSees.get_fiber_data import *

Frame_number=9
dict=Frame[str(Frame_number)]
Frame_details=Frame_data(dict)
Frame_details.geometric_imperfection_ratio=+1/500
if Frame_details.geometric_imperfection_ratio>0:
    wind_load_dirn='right'
else:
    wind_load_dirn='left'


Frame=Moment_Frame_2D(Frame_details.bay_width, Frame_details.story_height, Frame_details.column_no_of_ele, Frame_details.beam_no_of_ele,
                beam_section=Frame_details.beam_section,
                column_section=Frame_details.column_section,
                support=Frame_details.support,
                D_floor_intensity=Frame_details.D_floor_intensity,
                D_roof_intensity=Frame_details.D_roof_intensity,
                L_floor_intensity=Frame_details.L_floor_intensity,
                L_roof_intensity=Frame_details.L_roof_intensity,
                Wind_load_floor=Frame_details.Wind_load_floor,
                Wind_load_roof=Frame_details.Wind_load_roof,
                Wall_load=Frame_details.Wall_load,
                load_combination_multipliers=Frame_details.load_comb_multipliers,
                Frame_id=Frame_details.Frame_id,
                Residual_Stress=False,
                Inelastic_analysis=False,
                Second_order_effects=False,
                stiffness_reduction=0.8,
                nip=4,
                wind_load_dirn=wind_load_dirn)

Frame.generate_Nodes_and_Element_Connectivity()
# print(Frame.Main_Nodes)
# print(Frame.NODES_TO_FIX)
# print(Frame.column_intermediate_nodes)
print(Frame.column_connectivity)
# print(Frame.beam_intermediate_nodes)
print(Frame.beam_connectivity)
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
Frame.create_distorted_nodes_and_element_connectivity(Frame_details.geometric_imperfection_ratio)
# print(Frame.Main_Nodes)
Frame.build_ops_model()
 
# Frame.plot_model()
# print(Frame.roof_beams)
# print(Frame.column_connectivity)
Frame.run_gravity_analysis(steps=1000,plot_defo=True)
Frame.save_moments_by_member()
# x,y,A,m=get_fiber_data(str(1),plot_fibers=True)
# print(x)
# print(y)
# print(A)
# print(m)
# m_int = list(map(int,m))
# plt.scatter(x,y,A,m_int)
# plt.gca().axis('equal')
# plt.show()

# Frame.plot_all_fiber_section_in_the_model()
# print(Frame.beam_connectivity)
# print(Frame.sorted_element_connectivity)
# print(Frame.member_list)
Frame.plot_model()
# Frame.display_node_coords()