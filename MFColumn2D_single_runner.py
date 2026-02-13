from MFColumn2D import *
from Plots import line_plot
from Ziemian_database import Frame_Info, convert_dict_items_to_class_attributes,Analysis_Info
from libdenavit.OpenSees.get_fiber_data import *
import opsvis as opsv 

Frame_number='Trial_Col'
Analysis_type='GMNIA'
analysis_folder = os.path.join(Frame_number, Analysis_type)
Frame_dict=Frame_Info[str(Frame_number)]
Frame_details=convert_dict_items_to_class_attributes(Frame_dict)
if Frame_details.geometric_imperfection_ratio>0:
    wind_load_dirn='right'
else:
    wind_load_dirn='left'

Analysis_dict=Analysis_Info[str(Analysis_type)]
Analysis_details=convert_dict_items_to_class_attributes(Analysis_dict)

Frame=MFColumn_2D(Frame_details.bay_width, Frame_details.story_height, Frame_details.column_no_of_ele, Frame_details.beam_no_of_ele,
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
                    Residual_Stress=Analysis_details.Residual_Stress,
                    Elastic_analysis=Analysis_details.Elastic_analysis,
                    Second_order_effects=Analysis_details.Second_order_effects,
                    stiffness_reduction=Analysis_details.stiffness_reduction,
                    strength_reduction=Analysis_details.strength_reduction,
                    Notional_load=Analysis_details.Notional_load,
                    Geometric_Imperfection=Analysis_details.Geometric_Imperfection,
                    geometric_imperfection_ratio=Frame_details.geometric_imperfection_ratio,
                    nip=3,
                    mat_type='Steel01',
                    wind_load_dirn=wind_load_dirn,
                    Leaning_column=Frame_details.Leaning_column,
                    Leaning_column_offset=Frame_details.Leaning_column_offset,
                    Leaning_column_floor_load=Frame_details.Leaning_column_floor_load,
                    Leaning_column_roof_load=Frame_details.Leaning_column_roof_load)

# print(Frame.get_del2_over_del1())
Frame.generate_Nodes_and_Element_Connectivity()
# print(Frame.Main_Nodes)
# print(Frame.NODES_TO_FIX)
# print(Frame.column_intermediate_nodes)
# print("Column_Connectivity",Frame.column_connectivity)
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
# print(Frame.axis_i_roof_nodes(2))
# input()
# print(Frame.bay_i_internal_floor_nodes(2))
# print(Frame.bay_i_internal_roof_nodes(2))
# print(Frame.bay_i_floor_load(1))
# print(Frame.bay_i_roof_load(2))



# Frame.plot_model()
Frame.create_distorted_nodes_and_element_connectivity()
# print(Frame.all_nodes)
# print(Frame.Main_Nodes)
Frame.build_ops_model()
# disp=Frame.run_load_controlled_analysis(steps=100,plot_defo=False)

opsv.plot_model()
Frame.plot_model()
# Frame.build_ops_model()
# Frame.add_dead_live_wind_wall_loads()
# target_disp=-10 if disp<0 else 10
results,_ =Frame.run_displacement_controlled_analysis(plot_defo=True,analysis='non_proportional_limit_point',vertical_load_scale=1,lateral_load_scale=1,control_dir='L')
Frame.plot_model()
# --- Plot 1: λ vs displacement ---

filename = os.path.join(analysis_folder, f'load_ratio_vs_disp_{Frame_number}_{Analysis_type}_V.png')
line_plot(results.control_node_displacement, results.load_ratio,
        xlabel='Displacement at Control Node', ylabel='Load Ratio λ',
        title='Load Ratio vs Displacement', filename=filename)

# --- Plot 2: λ vs base shear ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_base_shear_{Frame_number}_{Analysis_type}_V.png')
line_plot(results.base_shear, results.load_ratio,
        xlabel='Base Shear', ylabel='Load Ratio λ',
        title='Load Ratio vs Base Shear', filename=filename)

# --- Plot 3: λ vs vertical reaction ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_vertical_reaction_{Frame_number}_{Analysis_type}_V.png')
line_plot(results.vertical_reaction, results.load_ratio,
        xlabel='Vertical Reaction', ylabel='Load Ratio λ',
        title='Load Ratio vs Vertical Reaction', filename=filename)

# --- Plot 4: λ vs max tensile strain ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_strain_{Frame_number}_{Analysis_type}_V.png')
line_plot(results.absolute_maximum_strain, results.load_ratio,
        xlabel='Maximum Tensile Strain', ylabel='Load Ratio λ',
        title='Load Ratio vs Tensile Strain', filename=filename)

# --- Plot 5: eigenvalue vs λ ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_eigenvalue_{Frame_number}_{Analysis_type}_V.png')
line_plot(results.load_ratio, results.lowest_eigenvalue,
        xlabel='Load Ratio λ', ylabel='Lowest Eigenvalue',
        title='Eigenvalue vs Load Ratio', filename=filename)

# --- Plot 6: λ vs P_M_M_interaction ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_P_M_M_interaction_{Frame_number}_{Analysis_type}_V.png')
line_plot(results.load_ratio, results.max_P_M_M_interaction,
        xlabel='Load Ratio λ', ylabel='max_P_M_M_interaction',
        title='P_M_M_interaction vs Load Ratio', filename=filename)
# Frame.save_moments_by_member()

# Frame.plot_all_fiber_section_in_the_model()
# print(Frame.beam_connectivity)
# print(Frame.sorted_element_connectivity)
# print(Frame.member_list)

# Frame.display_node_coords()