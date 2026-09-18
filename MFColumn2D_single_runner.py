from MFColumn2D import *
from Plots import line_plot
from Columns_config import Frame_Info
from Analysis_config import Analysis_Info
from Materials_config import Material_Info
from libdenavit.OpenSees.get_fiber_data import *
import opsvis as opsv 
from helpers import m, Steel_Material,convert_dict_items_to_class_attributes,check_and_create_new_entries_in_column_config_file,save_analysis_history_figures
from libdenavit.OpenSees import plotting

json_wind_dirn_path="wind_load_dirn_data.json"


column_section_name='W8X31'
slenderness_ratios=50
bending_axes='x'
story_heights=get_storey_height_list_for_a_section(column_section_name,bending_axes,slenderness_ratios)
print(story_heights)
# input()
# story_height=round(story_height,4)
no_of_stories=2
support='f'
Analysis_type='GMNIA'
Material='50_ksi'  # Options: '36_ksi', '50_ksi'

'''
This function writes a new entry in the Frame_info dictionary in column_config.py if the entry does not already exist.
 To check, change the parameters and see a new entry made inside the file. This entry is later pulled as frame info for making 
 a MFColulmn_2D object. Other parameters related to the frame details (like floor load, roof load etc.) that are not passed to 
 the function are assigned default values.
 '''

frame_name=check_and_create_new_entries_in_column_config_file(column_section_name=column_section_name,
                                                              story_height=story_heights,
                                                              no_of_stories=no_of_stories,
                                                              bending_axes=bending_axes,
                                                              Material_type=Material,
                                                              Leaning_column=False,
                                                              Floor_to_Roof_load_ratio=1,
                                                              Leaning_Column_load_ratio=0,
                                                              support=support,)



'''
This chunk reads the parameters of the column from the Columns_config.py. You can check by printing some attributes of Frame_details.
This will print the values that we have in the Frame_info dictionary
'''
frame_key=str( frame_name)
Frame_dict=Frame_Info[frame_key]
Frame_details=convert_dict_items_to_class_attributes(Frame_dict)

print(Frame_details.column_section)
print(Frame_details.load_comb_multipliers)
print(Frame_details.D_floor_intensity)



'''
This chunk reads the details of the analysis form Analysis_config.py. You can check by printing some attributes of Analysis_details.
This will print the values that we have in the Frame_info dictionary
'''
Analysis_dict=Analysis_Info[str(Analysis_type)]
Analysis_details=convert_dict_items_to_class_attributes(Analysis_dict)
print(Analysis_details.Residual_Stress)
print(Analysis_details.Second_order_effects)


'''
This chunk reads the details of the material form Materials_config.py.
'''
Material_dict=Material_Info[str(Material)]
Material_details=convert_dict_items_to_class_attributes(Material_dict)
Steel=Steel_Material(mat_tag=1,E=Material_details.E,Fy=Material_details.Fy)



'''
This chunk checks if the lateral/wind loading direction (left or right) for the given frame details has already been reported from an OpenSees
analysis,
if yes, it reads the wind load dirn
else, it assigns None to wind load dirn
'''
wind_data=load_wind_dirn_data(json_wind_dirn_path=json_wind_dirn_path)
wind_data=ensure_frame_entry_exists(frame_key=frame_key,data=wind_data,json_wind_dirn_path=json_wind_dirn_path)
wind_load_dirn=wind_data[frame_key]["wind_load_dirn"]
initial_out_of_straightness_dirn=wind_data[frame_key]["initial_out_of_straightness_dirn"]
print(wind_load_dirn)
print(initial_out_of_straightness_dirn)

input()

'''
Making a MFColumn_2D object from the Frame_details, Analysis_details, Material_details read from the corresponding files.
'''
Frame=MFColumn_2D(Frame_details.bay_width, Frame_details.story_height, Frame_details.column_no_of_ele, Frame_details.beam_no_of_ele,  
                    beam_section=Frame_details.beam_section,
                    column_section=Frame_details.column_section,
                    support=Frame_details.support,
                    D_floor_intensity=Frame_details.D_floor_intensity,
                    D_roof_intensity=Frame_details.D_roof_intensity,
                    L_floor_intensity=Frame_details.L_floor_intensity,
                    L_roof_intensity=Frame_details.L_roof_intensity,
                #   Wind_load_floor=Frame_details.Wind_load_floor,
                    Base_Wind_load=Frame_details.Base_Wind_load,
                    Wall_load=Frame_details.Wall_load,
                    load_combination_multipliers=Frame_details.load_comb_multipliers,
                    Frame_id=Frame_details.Frame_id,
                    Material_obj=Steel,
                    Steel_Grade=Frame_details.Material_type,
                    Residual_Stress=Analysis_details.Residual_Stress,
                    Elastic_analysis=Analysis_details.Elastic_analysis,
                    Second_order_effects=Analysis_details.Second_order_effects,
                    stiffness_reduction=Analysis_details.stiffness_reduction,
                    strength_reduction=Analysis_details.strength_reduction,
                    Notional_load=Analysis_details.Notional_load,
                    Geometric_Imperfection=Analysis_details.Geometric_Imperfection,
                    geometric_imperfection_ratio=Frame_details.geometric_imperfection_ratio,
                    initial_out_of_straightness_ratio=Frame_details.initial_out_of_straightness_ratio,
                    nip=5,
                    mat_type='Steel01',
                    wind_load_dirn=wind_load_dirn,
                    initial_out_of_straightness_dirn=initial_out_of_straightness_dirn,
                    Leaning_column=Frame_details.Leaning_column,
                    Leaning_column_offset=Frame_details.Leaning_column_offset,
                    Leaning_column_floor_load=Frame_details.Leaning_column_floor_load,
                    Leaning_column_roof_load=Frame_details.Leaning_column_roof_load,
                    floor_nodes_free=False,
                    wind_load_same_for_all_h=False)



if Frame.wind_load_dirn is None:
        calculated_wind_load_dirn=Frame.get_lateral_loading_direction()   
        wind_data[frame_key]["wind_load_dirn"]=calculated_wind_load_dirn
        wind_data[frame_key]["wind_load_dirn_source"] = "analysis"
        save_wind_dirn_data(data=wind_data,json_wind_dirn_path=json_wind_dirn_path)


print("Wind Load Direction")
print(Frame.wind_load_dirn)
print(Frame.column_only_model)
input("Wind Load direction figured out")

if Frame.initial_out_of_straightness_dirn is None:
        calculated_initial_out_of_straightness_dirn=Frame.get_initial_out_of_straightness_direction()   
        wind_data[frame_key]["initial_out_of_straightness_dirn"]=calculated_initial_out_of_straightness_dirn
        wind_data[frame_key]["initial_out_of_straightness_dirn_source"] = "analysis"
        save_wind_dirn_data(data=wind_data,json_wind_dirn_path=json_wind_dirn_path)


print("Initial out of straightness")
print(Frame.initial_out_of_straightness_dirn)
input("Initial Out of Straightness figured out")

Frame.generate_Nodes_and_Element_Connectivity()
Frame.create_distorted_nodes_and_element_connectivity()
Frame.build_ops_model()
opsv.plot_model()
opsv.plot_load()
Frame.plot_model()
print(Frame.column_connectivity)
print(Frame.column_intermediate_nodes)
print(Frame.member_list)
print(Frame.sorted_column_connectivity)
print(Frame.column_member_list)
print("All nodes")
print(Frame.all_nodes)
input()


vertical_load_scale=0.2
lateral_load_scale=1
tolerance=1e-6
iterations=10
control_dir='L'
proportional_or_not='non_proportional_limit_point' 
'''
This is to compare the results with a vertical pushover vs a lateral pushover
'''
'''
This line makes a new folder with a name 'frame_name' inside 'Column_Results'
'''
analysis_folder = os.path.join('Column_Results', frame_name, Analysis_type,control_dir)
os.makedirs(analysis_folder, exist_ok=True)


results,fail_during_LCA =Frame.run_displacement_controlled_analysis(target_disp=1,steps=1000,plot_defo=False,num_steps_LCA=50, 
                                                      analysis=proportional_or_not,
                                                      vertical_load_scale=vertical_load_scale,
                                                      lateral_load_scale=lateral_load_scale,
                                                      control_dir=control_dir,try_smaller_steps=False,
                                                      live_plot=True,tolerance=tolerance,iterations=iterations) 
input()

fig, ax = plt.subplots(figsize=(10, 6))
opsv.plot_model(
    node_labels=False,
    element_labels=False,
    node_supports=True,
    ax=ax
)

opsv.plot_load(
    node_supports=True,
    ax=ax
)



ax.set_aspect("equal")
plt.show(block=True)


save_analysis_history_figures(
    parent_folder=analysis_folder,
    results=results,
    frame_name=frame_name,
    analysis_type=Analysis_type,
    ALR_H=lateral_load_scale,
    ALR_V=vertical_load_scale,
    tolerance=tolerance,
    iterations=iterations,
    control_dir=control_dir,
)

if proportional_or_not=='proportional_limit_point':

    del2_over_del1=Frame.get_del2_over_del1(vertical_load_scale=results.maximum_load_ratio_at_limit_point*vertical_load_scale,
                                            lateral_load_scale=results.maximum_load_ratio_at_limit_point*lateral_load_scale+0.5)


elif proportional_or_not=='non_proportional_limit_point':
    del2_over_del1=Frame.get_del2_over_del1(vertical_load_scale=vertical_load_scale,
                                            lateral_load_scale=results.maximum_load_ratio_at_limit_point*lateral_load_scale+5)
      