from Moment_Frame_2D_Main import *
from Ziemian_database import Frame_Info,Analysis_Info
from Ziemian_database import convert_dict_items_to_class_attributes
from libdenavit.OpenSees.get_fiber_data import *
import gc
from Plotting import plot_single_bar,line_plot
import os
import seaborn as sns

def MF_2D_runner(Frame_number,Analysis_type,control_dir='L',lateral_load_scale=1,vertical_load_scale=1,ops_anlaysis='proportional_limit_point'):
    try:
        Frame_dict=Frame_Info[str(Frame_number)]
        Frame_details=convert_dict_items_to_class_attributes(Frame_dict)
        Frame_details.geometric_imperfection_ratio=+1/500
        if Frame_details.geometric_imperfection_ratio>0:
            wind_load_dirn='right'
        else:
            wind_load_dirn='left'

        Analysis_dict=Analysis_Info[str(Analysis_type)]
        Analysis_details=convert_dict_items_to_class_attributes(Analysis_dict)


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
                        Residual_Stress=Analysis_details.Residual_Stress,
                        Elastic_analysis=Analysis_details.Elastic_analysis,
                        Second_order_effects=Analysis_details.Second_order_effects,
                        stiffness_reduction=Analysis_details.stiffness_reduction,
                        strength_reduction=Analysis_details.strength_reduction,
                        Notional_load=Analysis_details.Notional_load,
                        Geometric_Imperfection=Analysis_details.Geometric_Imperfection,
                        nip=3,
                        mat_type='Steel01',
                        wind_load_dirn=wind_load_dirn)

        Frame.generate_Nodes_and_Element_Connectivity()
        Frame.create_distorted_nodes_and_element_connectivity(Frame_details.geometric_imperfection_ratio)
        Frame.build_ops_model()
        results,fail_during_LCA=Frame.run_displacement_controlled_analysis(target_disp=5,plot_defo=False,control_dir=control_dir,
                    lateral_load_scale=lateral_load_scale,vertical_load_scale=vertical_load_scale,analysis=ops_anlaysis)
        # Frame.plot_model()
        return results,Frame_details.Frame_id,fail_during_LCA

    finally:
        # --- CLEANUP ---
        ops.wipe()          # clear OpenSees model
        if 'Frame' in locals():
            Frame.__dict__.clear()
            del Frame          # delete Python object
        gc.collect()


Frame_number= [13]      # 'SP36H'  ,  'UP36H'  ,  'SP36L'  ,  'UP36L'
Analysis_type= ['Second_Order_Inelastic'  ,  'AISC_Direct_Modeling'  , 'AISC_Notional_Loads']
# Analysis_type= [  'AISC_Notional_Loads'  ]
control_dir=['V']  # 'L'  'V'

# for frame_number in Frame_number:

#     max_load_ratio=[]
#     analysis_type_labels=[]
#     for analysis_type in Analysis_type:
#         analysis_type_labels.append(analysis_type.replace("_", " "))
#         results,frame_id,fail_during_LCA=MF_2D_runner(Frame_number=frame_number,Analysis_type=analysis_type,control_dir='V',ops_anlaysis='proportional_limit_point')

#         #    ---- λ vs displacement
#         os.makedirs(frame_id, exist_ok=True)
#         filename=os.path.join(frame_id,f'load_ratio_vs_disp_{frame_number}_{analysis_type}_{control_dir}.png')
#         line_plot(results.control_node_displacement_absolute,results.load_ratio,
#             xlabel='Displacement at Control Node', ylabel='Load Ratio λ',title='Load Ratio vs.  Displacement',filename=filename)

#         #     --- λ vs base_shear
#         filename=os.path.join(frame_id,f'load_ratio_vs_base_shear_{frame_number}_{analysis_type}_{control_dir}.png')
#         line_plot(results.base_shear,results.load_ratio,
#             xlabel='Base_shear', ylabel='Load Ratio λ',title='Load Ratio vs. Base Shear',filename=filename)
        
#         #     --- λ vs vertical reaction
#         filename=os.path.join(frame_id,f'load_ratio_vs_Vertical_Reaction_{frame_number}_{analysis_type}_{control_dir}.png')
#         line_plot(results.vertical_reaction,results.load_ratio,
#             xlabel='Vertical Reaction', ylabel='Load Ratio λ',title='Load Ratio vs. Vertical Reaction',filename=filename)

#         # --- λ vs max tensile strain
#         filename=os.path.join(frame_id,f'load_ratio_vs_strain_{frame_number}_{analysis_type}_{control_dir}_.png')
#         line_plot(results.absolute_maximum_strain,results.load_ratio,
#             xlabel='Maximum Tensile Strain',ylabel='Load Ratio λ',title='Load Ratio vs. Tensile Strain',filename=filename)

#         # --- eigenvalue vs λ  (x = λ, y = lowest eigenvalue)
#         filename=os.path.join(frame_id,f'load_ratio_vs_eigenvalue_{frame_number}_{analysis_type}_{control_dir}.png')
#         line_plot(results.load_ratio, results.lowest_eigenvalue,
#             xlabel='Load Ratio λ',ylabel='Lowest Eigenvalue',title='Eigenvalue vs. Load Ratio',filename=filename)
        
#         filename=os.path.join(frame_id,f'load_ratio_vs_P_M_M_interaction_{frame_number}_{analysis_type}_{control_dir}.png')
#         line_plot(results.load_ratio,results.max_P_M_M_interaction,
#             xlabel='Load Ratio λ',ylabel='max_P_M_M_interaction',title='P_M_M_interaction vs. Load Ratio',filename=filename)
        
#         max_load_ratio.append(results.maximum_load_ratio_at_limit_point)

#     os.makedirs(frame_id, exist_ok=True)
#     filename=os.path.join(frame_id, f"{frame_number}_Barplot")
#     plot_single_bar(max_load_ratio,filename=filename, x_labels=analysis_type_labels,title=f'{frame_number}:Comparison of Load Ratio',ylabel='Load Ratio')

for frame_number in Frame_number:
    fig, ax = plt.subplots(figsize=(5,5))
    palette = sns.color_palette("pastel", len(Analysis_type))
    for j,analysis_type in enumerate(Analysis_type):
        ALR_H = []
        ALR_V = []

        results, frame_id, fail_during_LCA = MF_2D_runner(
            Frame_number=frame_number,
            Analysis_type=analysis_type,
            lateral_load_scale=0,
            control_dir='V',
            ops_anlaysis='proportional_limit_point'
        )

        ALR_V_max = results.maximum_load_ratio_at_limit_point
        os.makedirs(frame_id, exist_ok=True)
        filename = os.path.join(
            frame_id, f'ALR_H_vs_ALR_V_{frame_number}.png'
        )

        # Sweep load ratios
        for i in np.arange(0, 0.99, 0.1):
            print(f'Running vertical load scale {i}')
            results, frame_id, fail_during_LCA = MF_2D_runner(
                Frame_number=frame_number,
                Analysis_type=analysis_type,
                vertical_load_scale=i*ALR_V_max,
                control_dir='L',
                ops_anlaysis='non_proportional_limit_point'
            )
            if fail_during_LCA:
                ALR_H.append(0)
                ALR_V.append(ALR_V_max)
                break
            else:
                ALR_H.append(results.maximum_load_ratio_at_limit_point)
                ALR_V.append(i*ALR_V_max)

        for i in np.arange(0.99, 1.1, 0.01):
            print(f'Running vertical load scale {i}')
            results, frame_id, fail_during_LCA = MF_2D_runner(
                Frame_number=frame_number,
                Analysis_type=analysis_type,
                vertical_load_scale=i*ALR_V_max,
                control_dir='L',
                ops_anlaysis='non_proportional_limit_point'
            )
            if fail_during_LCA:
                ALR_H.append(0)
                ALR_V.append(ALR_V_max)
                break
            else:
                ALR_H.append(results.maximum_load_ratio_at_limit_point)
                ALR_V.append(i*ALR_V_max)
        # Plot each analysis type on the same axes
        line_plot(
            ALR_H, ALR_V,
            xlabel='ALR_H',
            ylabel='ALR_V',
            title=f'Frame {frame_number}: ALR_V vs ALR_H',
            filename=filename,
            ax=ax,
            label=analysis_type,
            linewidth=1,
            markersize=2,
            color=palette[j]
        )
