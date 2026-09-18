from matplotlib.lines import Line2D
from Moment_Frame_2D_Main import *
from MFColumn2D import *
from Columns_config import Frame_Info
from Analysis_config import Analysis_Info
from Materials_config import Material_Info
from libdenavit.OpenSees.get_fiber_data import *
from libdenavit.OpenSees import plotting
from libdenavit.interaction_diagram_2d import InteractionDiagram2d,cart2pol
from Plots import plot_single_bar,line_plot
import os
import seaborn as sns
from typing import List, Optional, Tuple
import re
from helpers import load_wind_dirn_data,save_wind_dirn_data,ensure_frame_entry_exists,convert_dict_items_to_class_attributes,Steel_Material,ft,kip,WF_Database
from helpers import get_storey_height_list_for_a_section,save_analysis_history_figures
import psutil
import time
import matplotlib.pyplot as plt
import json
from numbers import Real
from mpi4py import MPI
json_wind_dirn_path="wind_load_dirn_data.json"

python_process = psutil.Process(os.getpid())

def MF_2D_runner(Frame_number,Analysis_type,tolerance,iterations,control_dir='L',
                 lateral_load_scale=1,vertical_load_scale=1,ops_anlaysis='proportional_limit_point',step_size=1000):
    # try:
    frame_key=str(Frame_number)
    Frame_dict=Frame_Info[frame_key]
    Frame_details=convert_dict_items_to_class_attributes(Frame_dict)
    

    Analysis_dict=Analysis_Info[str(Analysis_type)]
    Analysis_details=convert_dict_items_to_class_attributes(Analysis_dict)

    # Material_dict=Material_Info[str(Material_type)]
    Material_dict=Material_Info[str(Frame_details.Material_type)]
    Material_details=convert_dict_items_to_class_attributes(Material_dict)
    Steel=Steel_Material(mat_tag=1,E=Material_details.E,Fy=Material_details.Fy)

    wind_data=load_wind_dirn_data(json_wind_dirn_path=json_wind_dirn_path)
    wind_data=ensure_frame_entry_exists(frame_key=frame_key,data=wind_data,json_wind_dirn_path=json_wind_dirn_path)
    wind_load_dirn=wind_data[frame_key]["wind_load_dirn"]
    initial_out_of_straightness_dirn=wind_data[frame_key]["initial_out_of_straightness_dirn"]


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
                        nip=3,
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

    if Frame.initial_out_of_straightness_dirn is None:
            calculated_initial_out_of_straightness_dirn=Frame.get_initial_out_of_straightness_direction()   
            wind_data[frame_key]["initial_out_of_straightness_dirn"]=calculated_initial_out_of_straightness_dirn
            wind_data[frame_key]["initial_out_of_straightness_dirn_source"] = "analysis"
            save_wind_dirn_data(data=wind_data,json_wind_dirn_path=json_wind_dirn_path)


    Frame.generate_Nodes_and_Element_Connectivity()
    Frame.create_distorted_nodes_and_element_connectivity()
    Frame.build_ops_model()
    # steps=step_size if lateral_load_scale==0 else 1000
    steps=step_size
    results,fail_during_LCA=Frame.run_displacement_controlled_analysis(target_disp=1,steps=steps,plot_defo=False,control_dir=control_dir,
                lateral_load_scale=lateral_load_scale,vertical_load_scale=vertical_load_scale,analysis=ops_anlaysis,tolerance=tolerance,iterations=iterations)
    
    # Get all attributes and their values as a dictionary
    attributes_dict = {
        attr: getattr(Frame, attr) 
        for attr in dir(Frame) 
        if not attr.startswith('__') and not callable(getattr(Frame, attr))  # Exclude special attributes and methods
    }

    return results,fail_during_LCA,Frame


def Bar_plot_comparison(Frame_number,Analysis_type,Material_type):
    for frame_number in Frame_number:
        max_load_ratio = []
        analysis_type_labels = []

        for analysis_type in Analysis_type:
            print(f'Running {analysis_type}')
            analysis_type_labels.append(analysis_type.replace("_", " "))

            # --- Run analysis ---
            results, fail_during_LCA,Frame = MF_2D_runner(
                Frame_number=frame_number,
                Analysis_type=analysis_type,
                control_dir='L',
                ops_anlaysis='proportional_limit_point'
            )
            # --- Create folder structure ---
            # frame_id is usually like "Frame_1"
            analysis_folder = os.path.join("Column_Results/Analysis_History", Frame.Frame_id, analysis_type)
            os.makedirs(analysis_folder, exist_ok=True)
            
            # ---make animation of PMM evolution throught the analysis
            filename = os.path.join(analysis_folder,f'PMM_Movie{frame_number}_{analysis_type}.gif')
            title=f'PM_History_{frame_number}_{analysis_type}'
            plotting.animate_PMM_evolution(results.P_M_M_interaction_all_elements,save_path=filename,fps=40,title=title  )


            # --- Plot 1: λ vs displacement ---
            results.control_node_displacement_absolute[:] = [abs(x) if x is not None else None
                                    for x in results.control_node_displacement]
            filename = os.path.join(analysis_folder, f'load_ratio_vs_disp_{frame_number}_{analysis_type}_V.png')
            line_plot(results.control_node_displacement_absolute, results.load_ratio,
                    xlabel='Displacement at Control Node', ylabel='Load Ratio λ',
                    title='Load Ratio vs Displacement', filename=filename)

            # --- Plot 2: λ vs base shear ---
            filename = os.path.join(analysis_folder, f'load_ratio_vs_base_shear_{frame_number}_{analysis_type}_V.png')
            line_plot(results.base_shear, results.load_ratio,
                    xlabel='Base Shear', ylabel='Load Ratio λ',
                    title='Load Ratio vs Base Shear', filename=filename)

            # --- Plot 3: λ vs vertical reaction ---
            filename = os.path.join(analysis_folder, f'load_ratio_vs_vertical_reaction_{frame_number}_{analysis_type}_V.png')
            line_plot(results.vertical_reaction, results.load_ratio,
                    xlabel='Vertical Reaction', ylabel='Load Ratio λ',
                    title='Load Ratio vs Vertical Reaction', filename=filename)

            # --- Plot 4: λ vs max tensile strain ---
            filename = os.path.join(analysis_folder, f'load_ratio_vs_strain_{frame_number}_{analysis_type}_V.png')
            line_plot(results.absolute_maximum_strain, results.load_ratio,
                    xlabel='Maximum Tensile Strain', ylabel='Load Ratio λ',
                    title='Load Ratio vs Tensile Strain', filename=filename)

            # --- Plot 5: eigenvalue vs λ ---
            filename = os.path.join(analysis_folder, f'load_ratio_vs_eigenvalue_{frame_number}_{analysis_type}_V.png')
            line_plot(results.load_ratio, results.lowest_eigenvalue,
                    xlabel='Load Ratio λ', ylabel='Lowest Eigenvalue',
                    title='Eigenvalue vs Load Ratio', filename=filename)

            # --- Plot 6: λ vs P_M_M_interaction ---
            filename = os.path.join(analysis_folder, f'load_ratio_vs_P_M_M_interaction_{frame_number}_{analysis_type}_V.png')
            line_plot(results.load_ratio, results.max_P_M_M_interaction,
                    xlabel='Load Ratio λ', ylabel='max_P_M_M_interaction',
                    title='P_M_M_interaction vs Load Ratio', filename=filename)

            # --- Record max load ratio ---
            max_load_ratio.append(results.maximum_load_ratio_at_limit_point)

        # --- Barplot comparing all analysis types for this frame ---
        os.makedirs(os.path.join("Column_Results/Analysis_History", Frame.Frame_id), exist_ok=True)
        barplot_filename = os.path.join("Column_Results/Analysis_History", Frame.Frame_id, f"{frame_number}_Barplot.png")
        plot_single_bar(max_load_ratio, filename=barplot_filename,
                        x_labels=analysis_type_labels,
                        title=f'{frame_number}: Comparison of Load Ratio',
                        ylabel='Load Ratio',title_fontsize=16, label_fontsize=14, tick_fontsize=14, legend_fontsize=12,value_fontsize=12)

def Interaction_Plots(Frame_number,Analysis_type,proportional=False,plot=False,
                      save_pmm_plots=False, pmm_plot_dpi=300,
                      step_size_V=1000,step_size_L=1000,
                      tolerance_V=1e-10,tolerance_L=1e-8,iterations=10):
    def save_pmm_plot(pmm_values, filename):
        if not save_pmm_plots:
            return
        fig_pmm, _ = plotting.plot_PMM_Interaction_values(
            pmm_values,
            show=False
        )
        try:
            fig_pmm.savefig(filename, dpi=pmm_plot_dpi)
        finally:
            plt.close(fig_pmm)


    exit_msg_color_map = {
        'Analysis Failed': "green",
        'Eigenvalue Limit Reached': "red",
        'Extreme Steel Fiber Strain Limit Reached': "orange",
        'P_M_M interaction Limit Reached': "purple",
        'Analysis Failed In Load Controlled Loading before entering Displacement controlled Loading': "blue",
        'Moving to Displacement Controlled Analysis': "cyan",
        'Analysis Failed in Displacement Controlled Loading in Non Proportional Analysis':"black"
    }

    ##non proportional
    if not proportional:
        code="Non_Proportional"
        for frame_number in Frame_number:
            palette = sns.color_palette("tab10", len(Analysis_type))
            linestyles = ['-', '-', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 1))]

            fig_alr, ax_alr = None, None

            if plot:
                fig_alr, ax_alr = plt.subplots(figsize=(5, 5))

            # store intersection points for optional later use
            intersections = []   # list of (x, y, color)

            for j, analysis_type in enumerate(Analysis_type):
                ALR_H, ALR_V = [], []
                exit_message=[]

                results, fail_during_LCA,Frame = MF_2D_runner(
                    Frame_number=frame_number,
                    Analysis_type=analysis_type,
                    lateral_load_scale=0,
                    control_dir='V',
                    ops_anlaysis='proportional_limit_point',
                    step_size=step_size_V,
                    tolerance=tolerance_V,
                    iterations=iterations
                )


                analysis_folder = os.path.join("Column_Results/Analysis_History", Frame.Frame_id, analysis_type)
                os.makedirs(analysis_folder, exist_ok=True)
            
                # --- Initialize ALR values ---
                ALR_V_max = results.maximum_load_ratio_at_limit_point
                V=ALR_V_max
                H=0
                ALR_V.append(V)
                ALR_H.append(H)
                exit_message.append(results.exit_message)
                print(ALR_V_max)

                check_load_ratio_problem(results,frame_number,analysis_type,V,H)

                video_name = f"deformed_shape_video_ALR_V_{ALR_V[-1]}_ALR_H_{ALR_H[-1]}.mp4"
                output_video_path = os.path.join(analysis_folder, video_name)

                # create_video_from_frames_and_clear_folder(
                #     temporary_folder="temporary_folder",
                #     output_video_path=output_video_path,
                #     fps=10
                # )
                save_analysis_history_figures(analysis_folder,results,Frame.Frame_id, analysis_type,ALR_H=round(H,4),ALR_V=round(V,4),
                                                step_size=step_size_V,
                                                tolerance=tolerance_V,
                                                iterations=iterations)
                # --- Plot PMM interaction for base case ---
                save_pmm_plot(
                    results.P_M_M_interaction_all_elements[-1],
                    os.path.join(
                        analysis_folder,
                        f"PMM_{analysis_type}_ALRH_{ALR_H[0]:.2f}_ALRV_{ALR_V[0]:.2f}.png"
                    )
                )

                # --- Sweep vertical loads (0–0.8) ---
                for i in np.arange(0, 0.4, 0.02):
                    print(f"Running vertical load scale {i:.2f}")
                    # input()
                    results, fail_during_LCA,Frame = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='non_proportional_limit_point',
                        tolerance=tolerance_L,
                        step_size=step_size_L,
                        iterations=iterations
                    )
                    if fail_during_LCA:
                        V=i*ALR_V_max*results.maximum_load_ratio_at_limit_point
                        H=0
                        ALR_H.insert(-1, H)
                        ALR_V.insert(-1, V)
                        exit_message.insert(-1,results.exit_message) 
                        
                    else:
                        V=i * ALR_V_max
                        H= results.maximum_load_ratio_at_limit_point
                        ALR_H.insert(-1, H)
                        ALR_V.insert(-1, V)
                        exit_message.insert(-1,results.exit_message)  


                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )
                    save_analysis_history_figures(analysis_folder,results,Frame.Frame_id, analysis_type,ALR_H=round(H,4),ALR_V=round(V,4),
                                                  tolerance=tolerance_L,
                                                  step_size=step_size_L,
                                                    iterations=iterations)

                # --- Sweep vertical loads (0.8–1.0) ---
                for i in np.arange(0.4, 0.5, 0.05):
                    print(f"Running vertical load scale {i:.2f}")
                    # input()
                    results, fail_during_LCA,Frame = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='non_proportional_limit_point',
                        tolerance=tolerance_L,
                        step_size=step_size_L,
                        iterations=iterations
                    )
                    if fail_during_LCA:
                        V=i*ALR_V_max*results.maximum_load_ratio_at_limit_point
                        H=0
                        ALR_H.insert(-1, H)
                        ALR_V.insert(-1, V)
                        exit_message.insert(-1,results.exit_message) 
                    else:
                        V=i * ALR_V_max
                        H= results.maximum_load_ratio_at_limit_point
                        ALR_H.insert(-1, H)
                        ALR_V.insert(-1, V)
                        exit_message.insert(-1,results.exit_message)

                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )
                    save_analysis_history_figures(analysis_folder,results,Frame.Frame_id, analysis_type,ALR_H=round(H,4),ALR_V=round(V,4),
                                                  tolerance=tolerance_L,
                                                  step_size=step_size_L,
                                                iterations=iterations)

                for i in np.arange(0.5, 1.0, 0.1):
                    print(f"Running vertical load scale {i:.2f}")
                    # input()
                    results,  fail_during_LCA,Frame = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='non_proportional_limit_point',
                        tolerance=tolerance_L,
                        step_size=step_size_L,
                        iterations=iterations
                    )
                    if fail_during_LCA:
                        V=i*ALR_V_max*results.maximum_load_ratio_at_limit_point
                        H=0
                        ALR_H.insert(-1, H)
                        ALR_V.insert(-1, V)
                        exit_message.insert(-1,results.exit_message) 
                    else:
                        V=i * ALR_V_max
                        H= results.maximum_load_ratio_at_limit_point
                        ALR_H.insert(-1, H)
                        ALR_V.insert(-1, V)
                        exit_message.insert(-1,results.exit_message)

                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )
                    save_analysis_history_figures(analysis_folder,results,Frame.Frame_id, analysis_type,ALR_H=round(H,4),ALR_V=round(V,4),
                                                  tolerance=tolerance_L,
                                                  step_size=step_size_L,
                                                iterations=iterations)
                
                # Convert to arrays for intersection math
                ALR_H_arr = np.array(ALR_H, dtype=float)
                ALR_V_arr = np.array(ALR_V, dtype=float)

                if plot:
                # --- Plot this analysis type on the SAME axes ---
                    line_plot(
                        ALR_H_arr, ALR_V_arr,
                        xlabel='ALR_H',
                        ylabel='ALR_V',
                        ax=ax_alr,
                        label=analysis_type,
                        linewidth=1.1,
                        markersize=0.1,
                        color=palette[j % len(palette)],
                        linestyle=linestyles[j % len(linestyles)],
                        show=False
                    )

                    for x, y, msg in zip(ALR_H_arr, ALR_V_arr, exit_message):
                        ax_alr.scatter(
                            x, y,
                            color=exit_msg_color_map.get(msg, "gray"),
                            s=15,
                            alpha=0.80,          
                            edgecolors='black',  
                            linewidths=0.25,
                            zorder=5
                        )

            # Save figure for this frame in non-proportional case
            if plot:
                
                (x_cross,y_cross) = intersection_with_ray_from_origin(x=ALR_H, y=ALR_V, theta_deg=45)

                ax_alr.scatter(
                    x_cross, y_cross,
                    s=15,
                    facecolors='none',
                    edgecolors=palette[j % len(palette)],
                    linewidths=0.8,
                    zorder=5
                )

                ax_alr.plot(
                    [0, x_cross], [0, y_cross],
                    color='black',
                    linestyle='--',
                    linewidth=0.8,
                    alpha=0.7,
                    label='ALR_H = ALR_V'
                )

                ax_alr.set_xlim(left=0)
                ax_alr.set_ylim(bottom=0)

                ax_alr.set_title(f' {frame_number} ({code})')

                curve_legend = ax_alr.legend(
                    loc='upper right',
                    fontsize=7,
                    title='Analysis Type',
                    title_fontsize=8
                )
                ax_alr.add_artist(curve_legend)

                termination_handles = [
                    Line2D(
                        [0], [0],
                        marker='o',
                        color='none',
                        markerfacecolor=color,
                        markeredgecolor='black',
                        markeredgewidth=0.25,
                        markersize=5,
                        alpha=0.60,
                        label=msg
                    )
                    for msg, color in exit_msg_color_map.items()
                ]

                ax_alr.legend(
                    handles=termination_handles,
                    loc='lower left',
                    fontsize=6,
                    title='Termination Message',
                    title_fontsize=7,
                    frameon=True
                )

                fig_alr.tight_layout()
                
                # Ensure directory exists before saving
                os.makedirs(os.path.join("Column_Results/Analysis_History", Frame.Frame_id), exist_ok=True)
                fig_alr.savefig(
                    os.path.join("Column_Results/Analysis_History", Frame.Frame_id, f'{code}_ALR_H_vs_ALR_V_Frame{frame_number}_{analysis_type}.png'),
                    dpi=600
                )
                plt.close(fig_alr)

    else:
    ### proportional
        code="Proportional"
        for frame_number in Frame_number:
            palette = sns.color_palette("tab10", len(Analysis_type))
            linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 1))]


            fig_alr, ax_alr = None, None       
            if plot:
                fig_alr, ax_alr = plt.subplots(figsize=(5, 5))
            # NEW: store intersections if you want them later (optional)
            intersections = []

            for j, analysis_type in enumerate(Analysis_type):
                ALR_H, ALR_V = [], []
                exit_message=[]

                # --- Base case: vertical-controlled analysis ---
                results,  fail_during_LCA,Frame = MF_2D_runner(
                    Frame_number=frame_number,
                    Analysis_type=analysis_type,
                    lateral_load_scale=0,
                    control_dir='V',
                    ops_anlaysis='proportional_limit_point',
                    step_size=step_size_V,
                    tolerance=tolerance_V,
                    iterations=iterations
                )

                analysis_folder = os.path.join("Column_Results/Analysis_History", Frame.Frame_id, analysis_type)
                os.makedirs(analysis_folder, exist_ok=True)

                

                # --- Initialize ALR values ---
                ALR_V_max = results.maximum_load_ratio_at_limit_point
                V=ALR_V_max
                H=0
                ALR_V.append(V)
                ALR_H.append(H)
                exit_message.append(results.exit_message)

                check_load_ratio_problem(results,frame_number,analysis_type,V,H)

                # --- Plot PMM interaction for base case ---
                save_pmm_plot(
                    results.P_M_M_interaction_all_elements[-1],
                    os.path.join(
                        analysis_folder,
                        f"PMM_{analysis_type}_ALRH_{ALR_H[0]:.2f}_ALRV_{ALR_V[0]:.2f}.png"
                    )
                )
                save_analysis_history_figures(analysis_folder,results,Frame.Frame_id, analysis_type,ALR_H=round(H,4),ALR_V=round(V,4),
                                                  tolerance=tolerance_V,
                                                  step_size=step_size_V,
                                                iterations=iterations)

                # --- Sweep vertical loads: fine steps near 0 ---
                for i in np.arange(0, 0.1, 0.05):
                    
                    print(f"Running vertical load scale {i:.3f}")
                    
                    results, fail_during_LCA,Frame = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point',
                        step_size=step_size_L,
                        tolerance=tolerance_L,
                        iterations=iterations
                    )

                    if fail_during_LCA and results.maximum_load_ratio_at_limit_point < 0.01:
                        break
                    else:
                        H=results.maximum_load_ratio_at_limit_point
                        V=i * ALR_V_max * results.maximum_load_ratio_at_limit_point
                        ALR_H.insert(-1, H)
                        # proportional case: both scaled by the same factor λ
                        ALR_V.insert(-1, V)
                        exit_message.insert(-1,results.exit_message)

                        check_load_ratio_problem(results,frame_number,analysis_type,V,H)
                        # --- Plot PMM every 2nd step ---
                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_proportional_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )
                    save_analysis_history_figures(analysis_folder,results,Frame.Frame_id, analysis_type,ALR_H=round(H,4),ALR_V=round(V,4),
                                                  tolerance=tolerance_L,
                                                  step_size=step_size_L,
                                                iterations=iterations)

                # --- Sweep vertical loads: 0.1 to 1.0 ---
                for i in np.arange(0.1, 1.0, 0.1):
                    
                    print(f"Running vertical load scale {i:.2f}")
                    
                    results, fail_during_LCA,Frame= MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point',
                        step_size=step_size_L,
                        tolerance=tolerance_L,
                        iterations=iterations
                    )

                    if fail_during_LCA or results.maximum_load_ratio_at_limit_point < 0.01:
                        break
                    else:
                        H=results.maximum_load_ratio_at_limit_point
                        V=i * ALR_V_max * results.maximum_load_ratio_at_limit_point
                        ALR_H.insert(-1, H)
                        # proportional case: both scaled by the same factor λ
                        ALR_V.insert(-1, V)
                        exit_message.insert(-1,results.exit_message)

                        check_load_ratio_problem(results,frame_number,analysis_type,V,H)
                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_proportional_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )
                    save_analysis_history_figures(analysis_folder,results,Frame.Frame_id, analysis_type,ALR_H=round(H,4),ALR_V=round(V,4),
                                                  tolerance=tolerance_L,
                                                  step_size=step_size_L,
                                                iterations=iterations)

                # --- Sweep vertical loads: 1.0 to 4.0 ---
                for i in np.arange(1.0, 2, 0.1):
                    
                    print(f'{analysis_type}')
                    print(f"Running vertical load scale {i:.2f}")
                    # input()

    
                    
                    results,  fail_during_LCA,Frame = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point',
                        step_size=step_size_L,
                        tolerance=tolerance_L,
                        iterations=iterations
                    )

                    if fail_during_LCA or results.maximum_load_ratio_at_limit_point < 0.01:
                        break
                    else:
                        H=results.maximum_load_ratio_at_limit_point
                        V=i * ALR_V_max * results.maximum_load_ratio_at_limit_point
                        ALR_H.insert(-1, H)
                        # proportional case: both scaled by the same factor λ
                        ALR_V.insert(-1, V)
                        exit_message.insert(-1,results.exit_message)

                        check_load_ratio_problem(results,frame_number,analysis_type,V,H)
                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_proportional_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )
                    save_analysis_history_figures(analysis_folder,results,Frame.Frame_id, analysis_type,ALR_H=round(H,4),ALR_V=round(V,4),
                                                  tolerance=tolerance_L,
                                                  step_size=step_size_L,
                                                iterations=iterations)

                # ============================================================
                #   PLOTTING ALR_H vs ALR_V FOR THIS ANALYSIS TYPE
                # ============================================================

                ALR_H_arr = np.array(ALR_H, dtype=float)
                ALR_V_arr = np.array(ALR_V, dtype=float)

                if plot:
                    line_plot(
                        ALR_H_arr, ALR_V_arr,
                        xlabel='ALR_H',
                        ylabel='ALR_V',
                        ax=ax_alr,
                        label=analysis_type,
                        linewidth=1.0,
                        markersize=0.3,
                        color=palette[j % len(palette)],
                        linestyle=linestyles[j % len(linestyles)],
                        show=False
                    )

                    for x, y, msg in zip(ALR_H_arr, ALR_V_arr, exit_message):
                        ax_alr.scatter(
                            x, y,
                            color=exit_msg_color_map.get(msg, "gray"),
                            s=28,
                            alpha=0.80,          # makes markers faint
                            edgecolors='black',  # helps visibility
                            linewidths=0.25,
                            zorder=5
                        )

            # Save figure for this frame in proportional case
            if plot:
                (x_cross,y_cross) = intersection_with_ray_from_origin(x=ALR_H, y=ALR_V, theta_deg=45)

                ax_alr.scatter(
                    x_cross, y_cross,
                    s=15,
                    facecolors='none',
                    edgecolors=palette[j % len(palette)],
                    linewidths=0.8,
                    zorder=5
                )

                ax_alr.plot(
                    [0, x_cross], [0, y_cross],
                    color='black',
                    linestyle='--',
                    linewidth=0.8,
                    alpha=0.7,
                    label='ALR_H = ALR_V'
                )

                ax_alr.set_xlim(left=0)
                ax_alr.set_ylim(bottom=0)


                ax_alr.set_title(f'Frame {frame_number}: ALR_V vs ALR_H ({code})')

                curve_legend = ax_alr.legend(
                    loc='upper right',
                    fontsize=7,
                    title='Analysis Type',
                    title_fontsize=8
                )
                ax_alr.add_artist(curve_legend)

                termination_handles = [
                    Line2D(
                        [0], [0],
                        marker='o',
                        color='none',
                        markerfacecolor=color,
                        markeredgecolor='black',
                        markeredgewidth=0.25,
                        markersize=5,
                        alpha=0.60,
                        label=msg
                    )
                    for msg, color in exit_msg_color_map.items()
                ]

                ax_alr.legend(
                    handles=termination_handles,
                    loc='lower left',
                    fontsize=6,
                    title='Termination Message',
                    title_fontsize=7,
                    frameon=True
                )

                fig_alr.tight_layout()
                
                # Ensure directory exists before saving
                os.makedirs(os.path.join("Column_Results/Analysis_History", Frame.Frame_id), exist_ok=True)
                fig_alr.savefig(
                    os.path.join("Column_Results/Analysis_History", Frame.Frame_id, f'{code}_ALR_H_vs_ALR_V_Frame{frame_number}_{analysis_type}.png'),
                    dpi=600
                )
                plt.close(fig_alr)

    # plotting.plot_sfd()
    # plotting.plot_bmd()
    # plotting.plot_afd(scale=0.001)
    return ALR_H,ALR_V,Frame

def intersection_with_ray_from_origin(
    x: List[float],
    y: List[float],
    theta_deg: float,
    eps: float = 1e-10
) -> Optional[Tuple[float, float]]:
    """
    Intersect polyline (x[i],y[i]) with the ray from origin at angle theta_deg.

    - theta=0°  : returns a point on the curve with y=0 (typically the x-axis endpoint).
    - theta=90° : returns a point on the curve with x=0 (typically the y-axis endpoint).
    - else      : intersects with y = tan(theta)*x, x>=0.

    Returns (xi, yi) or None if not found.
    """

    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")
    if len(x) < 2:
        return None
    if not (0.0 - eps <= theta_deg <= 90.0 + eps):
        raise ValueError("theta_deg must be between 0 and 90 degrees (inclusive).")

    # --- Handle theta = 0 (x-axis): look for y = 0 on the curve ---
    if abs(theta_deg - 0.0) <= eps:
        candidates = [(xi, yi) for xi, yi in zip(x, y) if abs(yi) <= eps and xi >= -eps]
        if not candidates:
            return None
        # pick the farthest on +x (usually the endpoint like (max_x, 0))
        return max(candidates, key=lambda p: p[0])

    # --- Handle theta = 90 (y-axis): look for x = 0 on the curve ---
    if abs(theta_deg - 90.0) <= eps:
        candidates = [(xi, yi) for xi, yi in zip(x, y) if abs(xi) <= eps and yi >= -eps]
        if not candidates:
            return None
        # pick the farthest on +y (usually the endpoint like (0, max_y))
        return max(candidates, key=lambda p: p[1])

    # --- General case: 0 < theta < 90 ---
    m = math.tan(math.radians(theta_deg))  # slope of the ray

    best_x = None
    best_pt = None

    for i in range(len(x) - 1):
        Ax, Ay = x[i], y[i]
        Bx, By = x[i + 1], y[i + 1]
        dx, dy = (Bx - Ax), (By - Ay)

        # Solve intersection with y = m x along segment:
        # (Ay - m Ax) + u[(By-Ay) - m(Bx-Ax)] = 0
        c0 = Ay - m * Ax
        c1 = dy - m * dx

        if abs(c1) <= eps:
            # segment is (almost) parallel to the ray in this equation
            # if also c0 ~ 0 -> colinear (infinite intersections); skip
            continue

        u = -c0 / c1
        if -eps <= u <= 1.0 + eps:
            u = min(1.0, max(0.0, u))
            xi = Ax + u * dx
            yi = Ay + u * dy

            # on forward ray (first quadrant)
            if xi >= -eps and yi >= -eps:
                # closest intersection to origin along the ray ~ minimize xi (since cos(theta)>0)
                if best_x is None or xi < best_x:
                    best_x = xi
                    best_pt = (xi, yi)

    return best_pt


def report_usage(iteration, log_file="memory_usage.txt"):
    python_gb = python_process.memory_info().rss / 1024**3

    vscode_bytes = 0

    for process in psutil.process_iter(["name", "memory_info"]):
        try:
            process_name = process.info["name"]

            if (
                process_name
                and process_name.lower() == "code.exe"
            ):
                vscode_bytes += process.info["memory_info"].rss

        except (
            psutil.NoSuchProcess,
            psutil.AccessDenied,
            psutil.ZombieProcess
        ):
            pass

    vscode_gb = vscode_bytes / 1024**3
    open_figures = len(plt.get_fignums())

    message = (
        f"Iteration {iteration}: "
        f"Python={python_gb:.2f} GB, "
        f"VS Code={vscode_gb:.2f} GB, "
        f"Open figures={open_figures}"
    )

    # Display in the terminal
    print(message)

    # Append to the text file
    with open(log_file, "a", encoding="utf-8") as file:
        file.write(message + "\n")

def extract_metadata(duplicate_frame):
    details = {
        attr: getattr(duplicate_frame, attr)
        for attr in ['_init_spec']
        if not attr.startswith('__') and not callable(getattr(duplicate_frame, attr))
    }
    init_spec = details['_init_spec']
    kwargs = init_spec['kwargs']


    return {
        "Frame_id": init_spec['Frame_id'],
        "width_of_bay": init_spec['width_of_bay'],
        "storey_height": init_spec['storey_height'],
        "no_of_elements_column": init_spec['no_of_elements_column'],
        "no_of_elements_beam": init_spec['no_of_elements_beam'],

        "beam_section": init_spec['beam_section'],
        "column_section": init_spec['column_section'],
        "column_section_name": init_spec['column_section']['common_and_exceptions']['common'][0],
        "bending_axes": init_spec['column_section']['common_and_exceptions']['common'][1],
        "no_of_stories": len(init_spec['storey_height']),

        "support": kwargs['support'],

        "D_floor_intensity": kwargs['D_floor_intensity'],
        "D_roof_intensity": kwargs['D_roof_intensity'],
        "L_floor_intensity": kwargs['L_floor_intensity'],
        "L_roof_intensity": kwargs['L_roof_intensity'],
        "Base_Wind_load": kwargs['Base_Wind_load'],
        "Wall_load": kwargs['Wall_load'],

        "load_combination_multipliers": init_spec['load_combination_multipliers'],

        # "Material_obj": init_spec['Material_obj'],
        "Steel_Grade": init_spec['Steel_Grade'],

        "Residual_Stress": kwargs['Residual_Stress'],
        "Elastic_analysis": kwargs['Elastic_analysis'],
        "Second_order_effects": kwargs['Second_order_effects'],
        "stiffness_reduction": kwargs['stiffness_reduction'],
        "strength_reduction": kwargs['strength_reduction'],
        "Notional_load": kwargs['Notional_load'],
        "Geometric_Imperfection": kwargs['Geometric_Imperfection'],
        "geometric_imperfection_ratio": kwargs['geometric_imperfection_ratio'],
        "initial_out_of_straightness_ratio":kwargs['initial_out_of_straightness_ratio'],
        "initial_out_of_straightness_dirn":kwargs['initial_out_of_straightness_dirn'],

        "nip": kwargs['nip'],
        "mat_type": kwargs['mat_type'],
        "wind_load_dirn": kwargs['wind_load_dirn'],

        "Leaning_column": kwargs['Leaning_column'],
        "Leaning_column_offset": kwargs['Leaning_column_offset'],
        "Leaning_column_floor_load": kwargs['Leaning_column_floor_load'],
        "Leaning_column_roof_load": kwargs['Leaning_column_roof_load'],

        "floor_nodes_free": kwargs['floor_nodes_free'],
        "wind_load_same_for_all_h": kwargs['wind_load_same_for_all_h']

    }

def compute_results(duplicate_frame,interaction_curve,analysis, theta_list, calculate_del2_over_del1=True):
    results={
             'Theta':[],
             'ALR_H':[],
             'ALR_V':[],
             'del2_over_del1':[]
             }
    results['Original_ALR_H']=interaction_curve.idx
    results['Original_ALR_V']=interaction_curve.idy

    if calculate_del2_over_del1:

        for theta in theta_list:
            pathX, pathY = ((0, 0), (0, 1)) if theta == 90 else ((0, 1), (0, math.tan(math.radians(theta))))
            pt = interaction_curve.find_intersection(pathX, pathY)
            vertical_load_scale = pt[1] if pt is not None else None
            if vertical_load_scale==0:
                ult_lat_load_for_V0=pt[0]
                no_of_digits = len(str(int(abs(ult_lat_load_for_V0)))) if ult_lat_load_for_V0!=0 else 1
            lateral_load_scale = ult_lat_load_for_V0*(10**(-no_of_digits - 2)) if pt[0]<10e-6 else pt[0] if pt is not None else None
            del2_over_del1=duplicate_frame.get_del2_over_del1(vertical_load_scale=vertical_load_scale, lateral_load_scale=lateral_load_scale)
            results['Theta'].append(theta)
            results['ALR_H'].append(pt[0] if pt is not None else None)
            results['ALR_V'].append(pt[1] if pt is not None else None)
            results['del2_over_del1'].append(del2_over_del1)

    return results

def write_interaction_results_to_json_file(
    frame_list,
    analysis_list,
    theta_list,
    iterations,
    calculate_del2_over_del1=True,
    proportional: bool = False,
    theta_round: int = 6,
    step_size_V=1000,
    step_size_L=100,
    tolerance_V=1e-8,
    tolerance_L=1e-8
    
    ) :
    """
    Make a new json file for each unique configuration and store the results
    for different analysis choices inside the json file. Also store the metadata
    that contains the detials of the configuration.
    """

    # Normalize theta list for stable matching
    theta_list = [round(float(t), theta_round) for t in theta_list] 

    for i,frame in enumerate(frame_list):
        frame = str(frame)
        print(f"\n=== Frame: {frame} ===")
        if i % 5 == 0:
            report_usage(i)
            
        curves = {}
        for analysis in analysis_list:
            print(f"Running {analysis} for {frame}")
            try:
                ALR_H, ALR_V,duplicate_frame = Interaction_Plots(
                    Frame_number=[frame],
                    Analysis_type=[analysis],
                    proportional=proportional,
                    plot=True,
                    save_pmm_plots=False,
                    step_size_V=step_size_V,
                    step_size_L=step_size_L,
                    tolerance_V=tolerance_V,
                    tolerance_L=tolerance_L,
                    iterations=iterations
                )
            finally:
                plt.close("all")

            curves[analysis] = (ALR_H, ALR_V)

        metadata = extract_metadata(duplicate_frame)
        new_data = {"data": metadata}

        file_name = f"{metadata['Frame_id']}.json"
        output_dir = "Column_Results/json_files"
        os.makedirs(output_dir, exist_ok=True)
        file_path = os.path.join(output_dir, file_name)

        update_or_create_json(file_path, new_data)


        for analysis in analysis_list:
            ALR_H, ALR_V = curves[analysis]
            interaction_curve=InteractionDiagram2d(ALR_H, ALR_V)
            # interaction_curve.plot()
            # plt.legend()
            # plt.title(f'Interaction Diagram for Frame {frame}')
            # plt.xlabel('ALR_H')
            # plt.ylabel('ALR_V')
            # plt.grid()
            interaction_diagram_save_path=os.path.join("Column_Results/Analysis_History", f"{frame}")
            os.makedirs(interaction_diagram_save_path, exist_ok=True)
            # plt.savefig(os.path.join(interaction_diagram_save_path, f'interaction_diagram_before_writing_to_file_{analysis}_{frame}.png'))
            # plt.close()
            
            results = compute_results(duplicate_frame,interaction_curve,analysis,theta_list,calculate_del2_over_del1=calculate_del2_over_del1)
            original_results={"Original_ALR_H":ALR_H,
                              "Original_ALR_V":ALR_V}
            Original={f"{analysis}": original_results}
            Results = {f"{analysis}": results}
            # update_or_create_json(file_path,Original)
            update_or_create_json(file_path, Results)




if __name__ == "__main__":


    comm=MPI.COMM_WORLD
    pid=comm.Get_rank()
    n_process=comm.Get_size()


    step_size_V=100000
    step_size_L=100000
    tolerance_V=1e-8
    tolerance_L=1e-3
    iterations=10
    new_analysis_run=True

    column_section_names=['W14X43','W14X120','W14X311','W14X730','W18X311','W21X62','W40X264','W40X392','W40X593']
    column_section_names=['W14X43','W14X120','W14X311','W14X730','W18X311','W21X62','W40X264','W40X392','W40X593']
    column_section_names=['W14X43','W14X120','W14X311','W14X730','W18X311','W21X62','W40X264','W40X392','W40X593']
    column_section_names=['W14X43']

    slenderness_ratios = np.arange(40, 55, 2)
    
    
    No_of_stories=[1]
    

    bending_axes=['x','y']
    # bending_axes=['y']

    Analysis_type= [  'GMNIA','GNA','GNA_Notional_Loads']  
    Analysis_type= [ 'GMNIA','GNA_Notional_Loads', 'GNA_no_stiffness_reduction']  
    Analysis_type= [ 'GMNIA','GNA','GNA_Notional_Loads']  
    # Analysis_type= [ 'GNA']
    Material_type=["50_ksi"]

    supports=['f','p']
    # supports=['p']

    theta_list = np.linspace(0, 90,91) 

task_count = 0
for support in supports:
    for material in Material_type:
        for bending_axis in bending_axes:
            for no_of_story in No_of_stories:
                for column_section_name in column_section_names:

                    story_heights = get_storey_height_list_for_a_section(
                        column_section_name,
                        bending_axis,
                        slenderness_ratios,
                    )

                    print(story_heights)



                    for story_height in story_heights:
                        story_height = round(float(story_height), 3)

                        # Distribute every complete combination among MPI ranks
                        if task_count % n_process == pid:
                            print(
                                f"Rank {pid}: task {task_count}, "
                                f"material={material}, "
                                f"axis={bending_axis}, "
                                f"stories={no_of_story}, "
                                f"section={column_section_name}, "
                                f"height={story_height}"
                            )

                            frame_name = (
                                check_and_create_new_entries_in_column_config_file(
                                    column_section_name,
                                    story_height,
                                    no_of_story,
                                    bending_axis,
                                    material,
                                    Leaning_column=False,
                                    Floor_to_Roof_load_ratio=1,
                                    Leaning_Column_load_ratio=0,
                                    support=support
                                )
                            )

                            write_interaction_results_to_json_file(
                                frame_list=[frame_name],
                                analysis_list=Analysis_type,
                                theta_list=theta_list,
                                proportional=False,
                                step_size_V=step_size_V,
                                step_size_L=step_size_L,
                                tolerance_V=tolerance_V,
                                tolerance_L=tolerance_L,
                                iterations=iterations,
                                calculate_del2_over_del1=False,
                            )

                        task_count += 1

