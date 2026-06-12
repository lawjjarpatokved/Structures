from MFColumn2D import *
from Plots import line_plot
from Columns_config import Frame_Info
from Analysis_config import Analysis_Info
from Materials_config import Material_Info
from libdenavit.OpenSees.get_fiber_data import *
import opsvis as opsv 
from helpers import m, Steel_Material,convert_dict_items_to_class_attributes,check_and_create_new_entries_in_column_config_file

json_wind_dirn_path="wind_load_dirn_data.json"
control_dir='V' 
'''
This is to compare the results with a vertical pushover vs a lateral pushover
'''

column_section_name='W8X31'
story_height=11*ft  ## Default units is m
no_of_stories=1
bending_axes='x'
Analysis_type='GNA'
Material='50_ksi'  # Options: '36_ksi', '50_ksi'

'''
This function writes a new entry in the Frame_info dictionary in column_config.py if the entry does not already exist.
 To check, change the parameters and see a new entry made inside the file. This entry is later pulled as frame info for making 
 a MFColulmn_2D object. Other parameters related to the frame details (like floor load, roof load etc.) that are not passed to 
 the function are assigned default values.
 '''

frame_name=check_and_create_new_entries_in_column_config_file(column_section_name=column_section_name,
                                                              story_height=story_height,
                                                              no_of_stories=no_of_stories,
                                                              bending_axes=bending_axes)

'''
This line makes a new folder with a name 'frame_name' inside 'Column_Results'
'''
analysis_folder = os.path.join('Column_Results', frame_name, Analysis_type,control_dir)
os.makedirs(analysis_folder, exist_ok=True)


'''
This chunk reads the parameters of the column from the Columns_config.py. You can check by printing some attributes of Frame_details.
This will print the values that we have in the Frame_info dictionary
'''
frame_key=str( frame_name)
Frame_dict=Frame_Info[frame_key]
Frame_details=convert_dict_items_to_class_attributes(Frame_dict)
print(Frame_details.column_section)
print(Frame_details.load_comb_multipliers)
# Frame_details.geometric_imperfection_ratio=0.002

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
wind_load_dirn='left'


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
                    Wind_load_floor=Frame_details.Wind_load_floor,
                    Wind_load_roof=Frame_details.Wind_load_roof,
                    Wall_load=Frame_details.Wall_load,
                    load_combination_multipliers=Frame_details.load_comb_multipliers,
                    Frame_id=Frame_details.Frame_id,
                    Material_obj=Steel,
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



if Frame.wind_load_dirn is None:
        calculated_wind_load_dirn=Frame.get_lateral_loading_direction()   
        wind_data[frame_key]["wind_load_dirn"]=calculated_wind_load_dirn
        wind_data[frame_key]["wind_load_dirn_source"] = "analysis"
        save_wind_dirn_data(data=wind_data,json_wind_dirn_path=json_wind_dirn_path)

print(Frame.wind_load_dirn)
input('Lateral loading figured out')
# a=Frame.get_del2_over_del1(lateral_load_scale=0.001,vertical_load_scale=0.87)
# # print(a)
# input()

# Frame=Frame2.rebuild_with_overrides(storey_height=[40 * ft, 28 * ft],                Residual_Stress=False,
#                 Elastic_analysis=True,
#                 Second_order_effects=True,
#                 stiffness_reduction=0.8,
#                 strength_reduction=1,
#                 Geometric_Imperfection=True,
#                 Notional_load=False)

# print(Frame.get_del2_over_del1())
Frame.generate_Nodes_and_Element_Connectivity()
Frame.create_distorted_nodes_and_element_connectivity()
Frame.build_ops_model()
opsv.plot_model()
Frame.plot_model()

results,_ =Frame.run_displacement_controlled_analysis(target_disp=1,steps=1000,plot_defo=True,analysis='proportional_limit_point',vertical_load_scale=1,lateral_load_scale=0,control_dir=control_dir) 
Frame.plot_model()
# --- Plot 1: λ vs displacement ---
print('load ratio')
print(results.load_ratio)
# print(results.control_node_displacement)
input()

filename = os.path.join(analysis_folder, f'load_ratio_vs_disp_{frame_name}_{Analysis_type}_{control_dir}.png')
line_plot(results.control_node_displacement, results.load_ratio,
        xlabel='Displacement at Control Node', ylabel='Load Ratio λ',
        title='Load Ratio vs Displacement', filename=filename)

# --- Plot 2: λ vs base shear ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_base_shear_{frame_name}_{Analysis_type}_{control_dir}.png')
line_plot(results.base_shear, results.load_ratio,
        xlabel='Base Shear', ylabel='Load Ratio λ',
        title='Load Ratio vs Base Shear', filename=filename)

# --- Plot 3: λ vs vertical reaction ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_vertical_reaction_{frame_name}_{Analysis_type}_{control_dir}.png')
line_plot(results.vertical_reaction, results.load_ratio,
        xlabel='Vertical Reaction', ylabel='Load Ratio λ',
        title='Load Ratio vs Vertical Reaction', filename=filename)

# --- Plot 4: λ vs max tensile strain ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_strain_{frame_name}_{Analysis_type}_{control_dir}.png')
line_plot(results.absolute_maximum_strain, results.load_ratio,
        xlabel='Maximum Tensile Strain', ylabel='Load Ratio λ',
        title='Load Ratio vs Tensile Strain', filename=filename)

# --- Plot 5: eigenvalue vs λ ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_eigenvalue_{frame_name}_{Analysis_type}_{control_dir}.png')
line_plot(results.load_ratio, results.lowest_eigenvalue,
        xlabel='Load Ratio λ', ylabel='Lowest Eigenvalue',
        title='Eigenvalue vs Load Ratio', filename=filename)

# --- Plot 6: λ vs P_M_M_interaction ---
filename = os.path.join(analysis_folder, f'load_ratio_vs_P_M_M_interaction_{frame_name}_{Analysis_type}_{control_dir}.png')
line_plot(results.load_ratio, results.max_P_M_M_interaction,
        xlabel='Load Ratio λ', ylabel='max_P_M_M_interaction',
        title='P_M_M_interaction vs Load Ratio', filename=filename)
