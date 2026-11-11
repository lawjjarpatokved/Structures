from MFColumn2D import MFColumn_2D # @todo avoid using * in imports
from Analysis_config import Analysis_Info
from Materials_config import Material_Info
from helpers import Steel_Material,save_analysis_history_figures
import os


Frame_id= "Trial_MFColumn_2D"
column_section_name='W8X31'
bending_axes='x'
story_height=3
no_of_stories=2
support='f'
Material='50_ksi'  # Options: '36_ksi', '50_ksi'
Analysis_type='GMNIA'
wind_load_dirn='left'
initial_out_of_straightness_dirn = 'left'


'''
This chunk reads the details of the analysis form Analysis_config.py. You can check by printing some attributes of Analysis_details.
This will print the values that we have in the Frame_info dictionary
'''
Analysis_dict=Analysis_Info[str(Analysis_type)]

'''
This chunk reads the details of the material form Materials_config.py.
'''
Material_dict=Material_Info[str(Material)]
Steel=Steel_Material(mat_tag=2,E=Material_dict['E'],Fy=Material_dict['Fy'])


details_dict = {
    "width_of_bay": [],
    "storey_height": [story_height]*no_of_stories,
    "no_of_elements_column": 4,
    "no_of_elements_beam": 4,  # @todo - needed for MFColumn2D?
    "beam_section": {
        "common_and_exceptions": {
            "common": "W27X84"
        }
    }, # @todo - needed for MFColumn2D?
    "column_section": {
        "common_and_exceptions": {
            "common": (column_section_name, bending_axes)
        }
    },
    "load_combination_multipliers": [1, 0, 0, 1],
    "Frame_id": Frame_id,
    "Material_obj": Steel,
    "Steel_Grade": "50_ksi",
    "support": "f",
    "D_floor_intensity": 1000.0,
    "D_roof_intensity": 1000.0,
    "L_floor_intensity": 0.0,
    "L_roof_intensity": 0.0,
    "Base_Wind_load": 1000.0,
    "Wall_load": 0, # @todo - maybe the default on this could be zero.
    "Residual_Stress": Analysis_dict["Residual_Stress"],
    "Elastic_analysis": Analysis_dict["Elastic_analysis"],
    "Second_order_effects": Analysis_dict["Second_order_effects"],
    "stiffness_reduction": Analysis_dict["stiffness_reduction"],
    "strength_reduction": Analysis_dict["strength_reduction"],
    "Notional_load": Analysis_dict["Notional_load"],
    "Geometric_Imperfection": Analysis_dict["Geometric_Imperfection"],
    "P_M_M_interaction_limit":Analysis_dict["P_M_M_interaction_limit"],
    "geometric_imperfection_ratio": 0.002,
    "initial_out_of_straightness_ratio": 0.001,
    "nip": 5,
    "mat_type": "Steel01",
    "wind_load_dirn": wind_load_dirn,
    "initial_out_of_straightness_dirn": initial_out_of_straightness_dirn,
    "Leaning_column": True,
    "Leaning_column_offset": 2,
    "Leaning_column_floor_load": 0.1,
    "Leaning_column_roof_load": 0.1,
    "floor_nodes_free": False,
    "wind_load_same_for_all_h": False,
    "column_only_model": True  # @todo - needed for MFColumn2D?
}

# Create MFColumn_2D Object
Frame=MFColumn_2D.from_dict(details_dict=details_dict)
Frame.generate_Nodes_and_Element_Connectivity()
Frame.create_distorted_nodes_and_element_connectivity()
Frame.build_ops_model()

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
analysis_folder = os.path.join('Column_Results', Frame_id, Analysis_type,control_dir)
os.makedirs(analysis_folder, exist_ok=True)
results,fail_during_LCA =Frame.run_displacement_controlled_analysis(target_disp=1,steps=1000,plot_defo=False,num_steps_LCA=50, 
                                                      analysis=proportional_or_not,
                                                      vertical_load_scale=vertical_load_scale,
                                                      lateral_load_scale=lateral_load_scale,
                                                      control_dir=control_dir,try_smaller_steps=False,
                                                      live_plot=False,tolerance=tolerance,iterations=iterations)  

save_analysis_history_figures(
    parent_folder=analysis_folder,
    results=results,
    frame_name=Frame_id,
    analysis_type=Analysis_type,
    ALR_H=lateral_load_scale,
    ALR_V=vertical_load_scale,
    tolerance=tolerance,
    iterations=iterations,
    control_dir=control_dir,
)

# if proportional_or_not=='proportional_limit_point':
#     del2_over_del1=Frame.get_del2_over_del1(vertical_load_scale=results.maximum_load_ratio_at_limit_point*vertical_load_scale,
#                                             lateral_load_scale=results.maximum_load_ratio_at_limit_point*lateral_load_scale+0.5)


# elif proportional_or_not=='non_proportional_limit_point':
#     del2_over_del1=Frame.get_del2_over_del1(vertical_load_scale=vertical_load_scale,
#                                             lateral_load_scale=results.maximum_load_ratio_at_limit_point*lateral_load_scale+5)
      