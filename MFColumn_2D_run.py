from matplotlib.pyplot import plot
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
json_wind_dirn_path="wind_load_dirn_data.json"


def MF_2D_runner(Frame_number,Analysis_type,Material_type,control_dir='L',lateral_load_scale=1,vertical_load_scale=1,ops_anlaysis='proportional_limit_point'):
    # try:
    frame_key=str(Frame_number)
    Frame_dict=Frame_Info[frame_key]
    Frame_details=convert_dict_items_to_class_attributes(Frame_dict)
    # if Frame_details.geometric_imperfection_ratio>0:
    #     wind_load_dirn='right'
    # else:
    #     wind_load_dirn='left'

    Analysis_dict=Analysis_Info[str(Analysis_type)]
    Analysis_details=convert_dict_items_to_class_attributes(Analysis_dict)

    Material_dict=Material_Info[str(Material_type)]
    Material_details=convert_dict_items_to_class_attributes(Material_dict)
    Steel=Steel_Material(mat_tag=1,E=Material_details.E,Fy=Material_details.Fy)

    wind_data=load_wind_dirn_data(json_wind_dirn_path=json_wind_dirn_path)
    wind_data=ensure_frame_entry_exists(frame_key=frame_key,data=wind_data,json_wind_dirn_path=json_wind_dirn_path)
    wind_load_dirn=wind_data[frame_key]["wind_load_dirn"]


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


    Frame.generate_Nodes_and_Element_Connectivity()
    Frame.create_distorted_nodes_and_element_connectivity()
    Frame.build_ops_model()
    steps=100000 if lateral_load_scale==0 else 1000
    results,fail_during_LCA=Frame.run_displacement_controlled_analysis(target_disp=1,steps=steps,plot_defo=False,control_dir=control_dir,
                lateral_load_scale=lateral_load_scale,vertical_load_scale=vertical_load_scale,analysis=ops_anlaysis)
    
    # Frame.plot_model()
    return results,Frame_details.Frame_id,fail_during_LCA,Frame

    # finally:
    #     # --- CLEANUP ---
    #     ops.wipe()          # clear OpenSees model
    #     if 'Frame' in locals():
    #         Frame.__dict__.clear()
    #         del Frame          # delete Python object
    #     gc.collect()

def Bar_plot_comparison(Frame_number,Analysis_type,Material_type):
    for frame_number in Frame_number:
        max_load_ratio = []
        analysis_type_labels = []

        for analysis_type in Analysis_type:
            print(f'Running {analysis_type}')
            analysis_type_labels.append(analysis_type.replace("_", " "))

            # --- Run analysis ---
            results, frame_id, fail_during_LCA,_ = MF_2D_runner(
                Frame_number=frame_number,
                Analysis_type=analysis_type,
                Material_type=Material_type,
                control_dir='L',
                ops_anlaysis='proportional_limit_point'
            )
            # --- Create folder structure ---
            # frame_id is usually like "Frame_1"
            analysis_folder = os.path.join("Column_Results", frame_id, analysis_type)
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
        os.makedirs(os.path.join("Column_Results", frame_id), exist_ok=True)
        barplot_filename = os.path.join("Column_Results", frame_id, f"{frame_number}_Barplot.png")
        plot_single_bar(max_load_ratio, filename=barplot_filename,
                        x_labels=analysis_type_labels,
                        title=f'{frame_number}: Comparison of Load Ratio',
                        ylabel='Load Ratio',title_fontsize=16, label_fontsize=14, tick_fontsize=14, legend_fontsize=12,value_fontsize=12)

def Interaction_Plots(Frame_number,Analysis_type,Material_type,proportional=False,plot=False,
                      save_pmm_plots=False, pmm_plot_dpi=300):
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

                results, frame_id, fail_during_LCA,Frame = MF_2D_runner(
                    Frame_number=frame_number,
                    Analysis_type=analysis_type,
                    Material_type=Material_type,
                    lateral_load_scale=0,
                    control_dir='V',
                    ops_anlaysis='proportional_limit_point'
                )

                analysis_folder = os.path.join("Column_Results", frame_id, analysis_type)
                os.makedirs(analysis_folder, exist_ok=True)

                # --- Initialize ALR values ---
                ALR_V_max = results.maximum_load_ratio_at_limit_point
                ALR_V.append(ALR_V_max)
                ALR_H.append(0)
                exit_message.append(results.exit_message)

                # --- Plot PMM interaction for base case ---
                save_pmm_plot(
                    results.P_M_M_interaction_all_elements[-1],
                    os.path.join(
                        analysis_folder,
                        f"PMM_{analysis_type}_ALRH_{ALR_H[0]:.2f}_ALRV_{ALR_V[0]:.2f}.png"
                    )
                )

                # --- Sweep vertical loads (0–0.8) ---
                for i in np.arange(0, 0.2, 0.05):
                    print(f"Running vertical load scale {i:.2f}")
                    # input()
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        Material_type=Material_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='non_proportional_limit_point'
                    )
                    if fail_during_LCA:
                        break
                    else:
                        ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                        ALR_V.insert(-1, i * ALR_V_max)
                        exit_message.append(results.exit_message)

                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )

                # --- Sweep vertical loads (0.8–1.0) ---
                for i in np.arange(0.2, 0.9, 0.1):
                    print(f"Running vertical load scale {i:.2f}")
                    # input()
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        Material_type=Material_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='non_proportional_limit_point'
                    )
                    if fail_during_LCA:
                        break
                    else:
                        ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                        ALR_V.insert(-1, i * ALR_V_max)
                        exit_message.append(results.exit_message)
                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )

                for i in np.arange(0.9, 1.0, 0.01):
                    print(f"Running vertical load scale {i:.2f}")
                    # input()
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        Material_type=Material_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='non_proportional_limit_point'
                    )
                    if fail_during_LCA:
                        break
                    else:
                        ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                        ALR_V.insert(-1, i * ALR_V_max)
                        exit_message.append(results.exit_message)
                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )

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
                    [0, 1], [0, 1],
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
                os.makedirs(os.path.join("Column_Results", frame_id), exist_ok=True)
                fig_alr.savefig(
                    os.path.join("Column_Results", frame_id, f'{code}_ALR_H_vs_ALR_V_Frame{frame_number}_{analysis_type}.png'),
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
                results, frame_id, fail_during_LCA,Frame = MF_2D_runner(
                    Frame_number=frame_number,
                    Analysis_type=analysis_type,
                    Material_type=Material_type,
                    lateral_load_scale=0,
                    control_dir='V',
                    ops_anlaysis='proportional_limit_point'
                )

                analysis_folder = os.path.join("Column_Results", frame_id, analysis_type)
                os.makedirs(analysis_folder, exist_ok=True)

                # --- Initialize ALR values ---
                ALR_V_max = results.maximum_load_ratio_at_limit_point
                ALR_V.append(ALR_V_max)
                ALR_H.append(0)
                exit_message.append(results.exit_message)

                # --- Plot PMM interaction for base case ---
                save_pmm_plot(
                    results.P_M_M_interaction_all_elements[-1],
                    os.path.join(
                        analysis_folder,
                        f"PMM_{analysis_type}_ALRH_{ALR_H[0]:.2f}_ALRV_{ALR_V[0]:.2f}.png"
                    )
                )

                # --- Sweep vertical loads: fine steps near 0 ---
                for i in np.arange(0, 0.1, 0.05):
                    
                    print(f"Running vertical load scale {i:.3f}")
                    
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        Material_type=Material_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point'
                    )

                    if fail_during_LCA and results.maximum_load_ratio_at_limit_point < 0.01:
                        break
                    else:
                        ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                        # proportional case: both scaled by the same factor λ
                        ALR_V.insert(-1, i * ALR_V_max * results.maximum_load_ratio_at_limit_point)
                        exit_message.append(results.exit_message)

                        # --- Plot PMM every 2nd step ---
                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_proportional_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )

                # --- Sweep vertical loads: 0.1 to 1.0 ---
                for i in np.arange(0.1, 1.0, 0.1):
                    
                    print(f"Running vertical load scale {i:.2f}")
                    
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        Material_type=Material_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point'
                    )

                    if fail_during_LCA or results.maximum_load_ratio_at_limit_point < 0.01:
                        break
                    else:
                        ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                        ALR_V.insert(-1, i * ALR_V_max * results.maximum_load_ratio_at_limit_point)
                        exit_message.append(results.exit_message)
                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_proportional_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )

                # --- Sweep vertical loads: 1.0 to 4.0 ---
                for i in np.arange(1.0, 40.0, 0.1):
                    
                    print(f'{analysis_type}')
                    print(f"Running vertical load scale {i:.2f}")
                    # input()

    
                    
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        Material_type=Material_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point'
                    )

                    if fail_during_LCA or results.maximum_load_ratio_at_limit_point < 0.01:
                        break
                    else:
                        ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                        ALR_V.insert(-1, i * ALR_V_max * results.maximum_load_ratio_at_limit_point)
                        exit_message.append(results.exit_message)
                        if int(i * 10) % 2 == 0:
                            save_pmm_plot(
                                results.P_M_M_interaction_all_elements[-1],
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_proportional_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                )
                            )

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
                    [0, 1], [0, 1],
                    color='black',
                    linestyle='--',
                    linewidth=0.8,
                    alpha=0.7,
                    label='ALR_H = ALR_V'
                )

                ax_alr.set_xlim(left=0)
                ax_alr.set_ylim(bottom=0)

                # ax_alr.set_title(f'Frame {frame_number}: ALR_V vs ALR_H ({code})')
                # ax_alr.legend()
                # fig_alr.tight_layout()

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
                os.makedirs(os.path.join("Column_Results", frame_id), exist_ok=True)
                fig_alr.savefig(
                    os.path.join("Column_Results", frame_id, f'{code}_ALR_H_vs_ALR_V_Frame{frame_number}_{analysis_type}.png'),
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

def write_interaction_results(
    csv_path: str,
    frame_list,
    analysis_list,
    Material_type,
    theta_list,
    proportional: bool = False,
    theta_round: int = 6
    ) -> pd.DataFrame:
    """
    Upsert interaction results into a CSV keyed by (Frame, Theta).

    Behavior:
    - If new theta values appear: new rows are added for those (Frame, Theta).
    - If new analyses appear: new columns are created (<analysis>_ALR_H, <analysis>_ALR_V).
    - Only columns for analyses in analysis_list are filled; other analysis columns remain NaN.
    """

    # Normalize theta list for stable matching
    theta_list = [round(float(t), theta_round) for t in theta_list]

    # --- Load existing or start new ---
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
    else:
        df = pd.DataFrame(columns=["Frame", "Theta"])

    # Ensure base columns exist
    if "Frame" not in df.columns:
        df["Frame"] = pd.Series(dtype="object")
    if "Theta" not in df.columns:
        df["Theta"] = pd.Series(dtype="float64")

    # Normalize existing theta + frame
    df["Frame"] = df["Frame"].astype(str)
    df["Theta"] = pd.to_numeric(df["Theta"], errors="coerce").round(theta_round)

    # --- Ensure required rows exist (Frame, Theta) ---
    existing_keys = set(zip(df["Frame"], df["Theta"]))

    new_rows = []
    for frame in frame_list:
        frame = str(frame)
        for theta in theta_list:
            if (frame, theta) not in existing_keys:
                new_rows.append({"Frame": frame, "Theta": theta})

    if new_rows:
        df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

    # --- Ensure required columns exist for the analyses you are running NOW ---
    for analysis in analysis_list:
        hcol = f"{analysis}_ALR_H"
        vcol = f"{analysis}_ALR_V"
        dcol = f"{analysis}_del2_over_del1"
        if hcol not in df.columns:
            df[hcol] = pd.Series(dtype='float64')
        if vcol not in df.columns:
            df[vcol] = pd.Series(dtype='float64')
        if dcol not in df.columns:
            df[dcol] = pd.Series(dtype='float64')

    # --- Index for updates ---
    df.set_index(["Frame", "Theta"], inplace=True)

    # --- Compute + fill only the analyses requested ---
    for frame in frame_list:
        frame = str(frame)
        print(f"\n=== Frame: {frame} ===")

        curves = {}
        for analysis in analysis_list:
            print(f"Running {analysis} for {frame}")
            # if analysis=='GMNIA':
                # input('Running GMNIA')
            try:
                ALR_H, ALR_V,duplicate_frame = Interaction_Plots(
                    Frame_number=[frame],
                    Analysis_type=[analysis],
                    Material_type=Material_type,
                    proportional=proportional,
                    plot=True,
                    save_pmm_plots=False
                )
            finally:
                plt.close("all")
            curves[analysis] = (ALR_H, ALR_V)
        print(theta_list)
        print(curves)
        # input()
        for analysis in analysis_list:
            ALR_H, ALR_V = curves[analysis]
            print(f"ALR_H: {ALR_H}")
            print(f"ALR_V: {ALR_V}")
            print(f"Plotting interaction diagram for {frame} - {analysis}")
            interaction_curve=InteractionDiagram2d(ALR_H, ALR_V)
            interaction_curve.plot()
            plt.legend()
            plt.title(f'Interaction Diagram for Frame {frame}')
            plt.xlabel('ALR_H')
            plt.ylabel('ALR_V')
            plt.grid()
            interaction_diagram_save_path=os.path.join("Column_Results", f"{frame}")
            os.makedirs(interaction_diagram_save_path, exist_ok=True)
            plt.savefig(os.path.join(interaction_diagram_save_path, f'interaction_diagram_before_writing_to_file_{analysis}_{frame}.png'))
            plt.close()
            # input()

            
            for theta in theta_list:
                
                print(theta)

                
                if theta==90:
                    pathX=(0,0)
                    pathY=(0,1)
                else:
                    pathX=(0,1)
                    pathY=(0,math.tan(math.radians(theta)))
                print(pathY)
                pt=interaction_curve.find_intersection(pathX,pathY)
                print(pt)
                # input()

                # pt = intersection_with_ray_from_origin(ALR_H, ALR_V, theta)
                # print(theta)
                # print(pt)
                # input()
                vertical_load_scale = pt[1] if pt is not None else None
                if vertical_load_scale==0:
                    ult_lat_load_for_V0=pt[0]
                    no_of_digits = len(str(int(abs(ult_lat_load_for_V0)))) if ult_lat_load_for_V0!=0 else 1
                    print(vertical_load_scale)
                    print(ult_lat_load_for_V0)
                    # input()
                print(ult_lat_load_for_V0)
                print(no_of_digits)
                print(10**(-no_of_digits - 2))
                # input(  )
                lateral_load_scale = ult_lat_load_for_V0*(10**(-no_of_digits - 2)) if pt[0]<10e-6 else pt[0] if pt is not None else None
                # lateral_load_scale = pt[0]
                del2_over_del1=duplicate_frame.get_del2_over_del1(vertical_load_scale=vertical_load_scale, lateral_load_scale=lateral_load_scale)

                # del2_over_del1=1

                hcol = f"{analysis}_ALR_H"
                vcol = f"{analysis}_ALR_V"
                dcol=f"{analysis}_del2_over_del1"

                if pt is None:
                    df.loc[(frame, theta), hcol] = np.nan
                    df.loc[(frame, theta), vcol] = np.nan
                    df.loc[(frame, theta), dcol] = np.nan
                else:
                    df.loc[(frame, theta), hcol] = float(pt[0])
                    df.loc[(frame, theta), vcol] = float(pt[1])
                    df.loc[(frame, theta), dcol] = float(del2_over_del1)

    # --- Save ---
    df.reset_index(inplace=True)
    df.sort_values(["Frame", "Theta"], inplace=True)
    
    # Convert all numeric columns to float to avoid text formatting warnings in Excel
    for col in df.columns:
        if col not in ['Frame', 'Theta']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    df.to_csv(csv_path, index=False)

    return df

def return_radial_error_betn_analyses(df: pd.DataFrame, analysis1: str, analysis2: str, tolerance=0.05) -> pd.Series:
    """
    Compare analysis1 against a tolerance-expanded analysis2.

    Positive value means analysis1 is within the expanded analysis2 overall.
    Negative value means analysis1 exceeds the expanded analysis2 overall.

    Here, analysis2 is enlarged by (1 + tolerance).
    """
    ana1_ALR_H = f"{analysis1}_ALR_H"
    ana1_ALR_V = f"{analysis1}_ALR_V"
    ana2_ALR_H = f"{analysis2}_ALR_H"
    ana2_ALR_V = f"{analysis2}_ALR_V"

    if (
        ana1_ALR_H not in df.columns or
        ana2_ALR_H not in df.columns or
        ana1_ALR_V not in df.columns or
        ana2_ALR_V not in df.columns
    ):
        raise ValueError(
            f"Columns '{ana1_ALR_H}', '{ana2_ALR_H}', '{ana1_ALR_V}', and/or '{ana2_ALR_V}' not found in DataFrame."
        )

    val1_H = pd.to_numeric(df[ana1_ALR_H], errors="coerce")
    val1_V = pd.to_numeric(df[ana1_ALR_V], errors="coerce")
    val2_H = pd.to_numeric(df[ana2_ALR_H], errors="coerce")
    val2_V = pd.to_numeric(df[ana2_ALR_V], errors="coerce")

    val2_H_tol = (1 + tolerance) * val2_H
    val2_V_tol = (1 + tolerance) * val2_V

    delta_H = val2_H_tol - val1_H
    delta_V = val2_V_tol - val1_V

    radial_error = delta_H + delta_V
    return radial_error

def plot_theta_vs_del2_over_del1(
    csv_path: str,
    frame_list: list,
    analyses_to_plot: list,
    theta_col: str = "Theta",
    frame_col: str = "Frame",
    show: bool = True,
    save_path: str | None = None
):
    """
    Plot Theta vs del2_over_del1 for selected analyses
    across multiple frames on one figure.
    """

    df = pd.read_csv(csv_path)

    plt.figure()

    for frame_name in frame_list:

        df_frame = df[df[frame_col].astype(str) == str(frame_name)].copy()

        if df_frame.empty:
            print(f"Warning: No data for Frame='{frame_name}'")
            continue

        df_frame[theta_col] = pd.to_numeric(df_frame[theta_col], errors="coerce")
        df_frame.sort_values(theta_col, inplace=True)

        for analysis in analyses_to_plot:
            col = f"{analysis}_del2_over_del1"

            if col not in df_frame.columns:
                print(f"Warning: Column '{col}' not found. Skipping.")
                continue

            y = pd.to_numeric(df_frame[col], errors="coerce")

            plt.plot(
                df_frame[theta_col],
                y,
                marker=".",
                linewidth=1,
                label=f"{frame_name} – {analysis}"
            )

    plt.xlabel("Theta (deg)")
    plt.ylabel("Δ₂ / Δ₁")
    plt.title("Theta vs Δ₂/Δ₁")
    plt.grid(True, alpha=0.3)
    plt.legend()

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    plt.close()

def plot_theta_vs_Radial_Errors(
    csv_path: str,
    frame_list: list,
    errors_to_plot: list,
    theta_col: str = "Theta",
    frame_col: str = "Frame",
    show: bool = True,
    save_path: str | None = None
):
    """
    Plot Theta vs Radial Error for selected error comparisons
    across multiple frames on one figure.
    """

    df = pd.read_csv(csv_path)

    # Compute error columns from the error names provided
    for error_name in errors_to_plot:
        print(error_name)
        start = error_name.find('(') + 1
        mid = error_name.find('<')
        first = error_name[start:mid]
    
        # Find text between < and )
        end = error_name.find(')')
        second = error_name[mid+1:end]


        col_name = f"Error({first}<{second})"
        df[col_name] = return_radial_error_betn_analyses(df, first, second)

    plt.figure()

    for frame_name in frame_list:

        df_frame = df[df[frame_col].astype(str) == str(frame_name)].copy()

        if df_frame.empty:
            print(f"Warning: No data for Frame='{frame_name}'")
            continue

        df_frame[theta_col] = pd.to_numeric(df_frame[theta_col], errors="coerce")
        df_frame.sort_values(theta_col, inplace=True)

        for error_name in errors_to_plot:
            start = error_name.find('(') + 1
            mid = error_name.find('<')
            first = error_name[start:mid]
            end = error_name.find(')')
            second = error_name[mid+1:end]
            
            col = f"Error({first}<{second})"

            if col not in df_frame.columns:
                print(f"Warning: Column '{col}' not found. Skipping.")
                continue

            y = pd.to_numeric(df_frame[col], errors="coerce")

            plt.plot(
                df_frame[theta_col],
                y,
                marker=".",
                linewidth=1,
                label=f"{frame_name} – {error_name}"
            )

    plt.xlabel("Theta (deg)")
    plt.ylabel("Radial Error")
    plt.title("Theta vs Radial Error")
    plt.grid(True, alpha=0.3)
    plt.legend()

    if save_path is not None:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    if show:
        plt.show()
    plt.close()

def parametric_h_min_radial_error(csv_path, column_section_name, bending_axes, no_of_stories, config_col='Frame', plot_slenderness=True):
    """
    Plot h vs min radial error or h/r vs min radial error.
    
    Parameters:
    -----------
    csv_path : str
        Path to the CSV file with results
    column_section_name : str
        Column section name (e.g., 'W27X84')
    bending_axes : str
        'x' or 'y' - which axis is being analyzed
    no_of_stories : int
        Number of stories
    config_col : str
        Column name in DataFrame to match
    plot_slenderness : bool
        If True, plot h/r (slenderness ratio). If False, plot h.
    """
    
    df = pd.read_csv(csv_path)
    plt.figure()
    ## Compute radial errors for all desired comparisons
    Radial_errors=[
           
            'Error(GNA_Notional_Loads<GNA)',
            'Error(GNA<GMNIA)',
            'Error(GNA_Notional_Loads<GMNIA)',
            
            ]
    # Compute error columns from the error names provided
    for error_name in Radial_errors:
        print(error_name)
        start = error_name.find('(') + 1
        mid = error_name.find('<')
        first = error_name[start:mid]
    
        # Find text between < and )
        end = error_name.find(')')
        second = error_name[mid+1:end]


        col_name = f"Error({first}<{second})"
        df[col_name] = return_radial_error_betn_analyses(df, first, second)

        # build prefix like W27X84_x_1X fromm the provided inputs
        prefix = f"{column_section_name}_{bending_axes}_{no_of_stories}X"    

        # --- filter matching rows ---
        df_sub = df[df[config_col].astype(str).str.startswith(prefix)].copy()

        if df_sub.empty:
            raise ValueError(f"No rows found matching prefix: {prefix}")

        # --- extract h from names like W27X84_x_1X18 ---
        df_sub["h"] = (
            df_sub[config_col]
            .astype(str)
            .str[len(prefix):]
        )
        df_sub["h"] = pd.to_numeric(df_sub["h"], errors="coerce")

        df_sub = df_sub.dropna(subset=["h"])

        if df_sub.empty:
            raise ValueError(f"Could not extract h values from {config_col} using prefix: {prefix}")

        # --- Get radius of gyration if plotting slenderness ratio ---
        if plot_slenderness:
            try:
                wf_db = WF_Database(column_section_name)
                r = wf_db.rx if bending_axes.lower() == 'x' else wf_db.ry
                df_sub["h_over_r"] = df_sub["h"] / r
                x_col = "h_over_r"
                x_label = f"h/r (Slenderness Ratio) [{column_section_name}, bending axis {bending_axes}]"
                filename_suffix = "h_over_r"
            except Exception as e:
                print(f"Warning: Could not fetch section properties: {e}")
                print(f"Falling back to h instead of h/r")
                x_col = "h"
                x_label = "h"
                filename_suffix = "h"
        else:
            x_col = "h"
            x_label = "h"
            filename_suffix = "h"

        

        for error_name in Radial_errors:
            if error_name not in df_sub.columns:
                print(f"Warning: {error_name} not found. Skipping.")
                continue

            y = pd.to_numeric(df_sub[error_name], errors="coerce")

            # ignore inf when taking minimum
            y = y.replace([np.inf, -np.inf], np.nan)

            temp = pd.DataFrame({
                x_col: df_sub[x_col],
                "error": y
            }).dropna()

            if temp.empty:
                print(f"Warning: No finite values found for {error_name}")
                continue

            grouped = temp.groupby(x_col, as_index=False)["error"].min()
            grouped = grouped.sort_values(x_col)

            plt.plot(
                grouped[x_col],
                grouped["error"],
                marker=".",
                linewidth=0.7,
                label=error_name
            )

        plt.xlabel(x_label)
        plt.ylabel("Minimum radial error")
        plt.title(f"{x_label} vs Minimum Radial Error ({prefix}*)")
        plt.grid(True, alpha=0.3)
        plt.legend()

        folder = os.path.join("Column_Results", "parametric h vs min radial error")
        os.makedirs(folder, exist_ok=True)

        save_path = os.path.join(folder, f"{prefix}_{filename_suffix}_vs_min_radial_error.png")
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

def parametric_slenderness_ratio_vs_min_radial_error(
    csv_path,
    config_col="Frame",
    save_folder=os.path.join("Column_Results", "parametric slenderness ratio vs min radial error")
):
    """
    Plot slenderness ratio (h/r) vs minimum radial error for all entries in the CSV.

    Expected frame naming format:
        SECTION_AXIS_STORIESXH
    Example:
        W27X84_x_1X18

    Meaning:
        section = W27X84
        axis = x
        no_of_stories = 1
        h = 18

    Parameters
    ----------
    csv_path : str
        Path to CSV file containing the results.
    config_col : str, optional
        Column name containing frame identifiers, by default "Frame".
    save_folder : str, optional
        Folder where the plot will be saved.

    Returns
    -------
    pd.DataFrame
        DataFrame with parsed frame information and computed h/r.
    """

    df = pd.read_csv(csv_path)

    radial_errors = [
        "Error(GNA_Notional_Loads<GNA)",
        "Error(GNA<GMNIA)",
        "Error(GNA_Notional_Loads<GMNIA)",
    ]

    # --------------------------------------------------
    # 1. Compute radial error columns if needed
    # --------------------------------------------------
    for error_name in radial_errors:
        start = error_name.find("(") + 1
        mid = error_name.find("<")
        end = error_name.find(")")

        first = error_name[start:mid]
        second = error_name[mid + 1:end]

        df[error_name] = return_radial_error_betn_analyses(df, first, second)

    # --------------------------------------------------
    # 2. Parse frame names
    # --------------------------------------------------
    # Expected pattern: SECTION_AXIS_STORIESXH
    # Example: W27X84_x_1X18
    pattern = re.compile(r"^(?P<section>.+?)_(?P<axis>[xyXY])_(?P<stories>\d+)X(?P<h>[-+]?\d*\.?\d+)$")

    parsed_rows = []

    for idx, frame_name in df[config_col].astype(str).items():
        match = pattern.match(frame_name.strip())
        if not match:
            print(f"Warning: Could not parse frame name: {frame_name}")
            continue

        section = match.group("section")
        axis = match.group("axis").lower()
        stories = int(match.group("stories"))
        h = float(match.group("h"))

        # Get radius of gyration
        try:
            wf_db = WF_Database(section)
            r = wf_db.rx if axis == "x" else wf_db.ry
        except Exception as e:
            print(f"Warning: Could not fetch WF properties for section {section}: {e}")
            continue

        if r is None or r == 0:
            print(f"Warning: Invalid radius of gyration for section {section}")
            continue

        row_dict = df.loc[idx].to_dict()
        row_dict["section"] = section
        row_dict["axis"] = axis
        row_dict["no_of_stories"] = stories
        row_dict["h"] = h
        row_dict["r"] = r
        row_dict["h_over_r"] = h / r

        parsed_rows.append(row_dict)

    if not parsed_rows:
        raise ValueError("No valid rows found after parsing frame names and section properties.")

    df_parsed = pd.DataFrame(parsed_rows)

    # --------------------------------------------------
    # 3. Plot min radial error vs h/r
    # --------------------------------------------------
    plt.figure()

    for error_name in radial_errors:
        if error_name not in df_parsed.columns:
            print(f"Warning: {error_name} not found in parsed DataFrame. Skipping.")
            continue

        y = pd.to_numeric(df_parsed[error_name], errors="coerce")
        y = y.replace([np.inf, -np.inf], np.nan)

        temp = pd.DataFrame({
            "h_over_r": pd.to_numeric(df_parsed["h_over_r"], errors="coerce"),
            "error": y
        }).dropna()

        if temp.empty:
            print(f"Warning: No finite values found for {error_name}")
            continue

        grouped = temp.groupby("h_over_r", as_index=False)["error"].min()
        grouped = grouped.sort_values("h_over_r")

        # plt.plot(
        #     grouped["h_over_r"],
        #     grouped["error"],
        #     marker=".",
        #     linewidth=0.7,
        #     label=error_name
        # )

        plt.scatter(
            grouped["h_over_r"],
            grouped["error"],
            label=error_name
        )

    plt.xlabel("h/r (Slenderness Ratio)")
    plt.ylabel("Minimum radial error")
    plt.title("Slenderness Ratio vs Minimum Radial Error")
    plt.grid(True, alpha=0.3)
    plt.legend()

    os.makedirs(save_folder, exist_ok=True)
    save_path = os.path.join(save_folder, "slenderness_ratio_vs_min_radial_error.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()

    return df_parsed

def add_structural_parameters_to_results(csv_path,config_col='Frame'):

    df = pd.read_csv(csv_path)



if __name__ == "__main__":
    new_analysis_run=True
    column_section_names=['W8X31','W14X43','W14X120','W14X311','W14X730','W18X311','W21X62','W40X264','W40X392','W40X593']
    column_section_names=['W14X43','W14X120','W14X311']
    story_heights=[7*ft,7.5*ft,8*ft,8.5*ft,9*ft,9.5*ft,10*ft,10.5*ft,11*ft,11.5*ft,12*ft,12.5*ft,13*ft,13.5*ft,14*ft,14.5*ft,15*ft,15.5*ft,16*ft,16.5*ft,17*ft,17.5*ft,18*ft,18.5*ft,19*ft,19.5*ft,20*ft]
    story_heights=[11*ft]
    story_heights=[7*ft,7.5*ft,8*ft,8.5*ft,9*ft,9.5*ft,10*ft,10.5*ft,11*ft,11.5*ft,12*ft,12.5*ft,13*ft,13.5*ft,14*ft]
    No_of_stories=[1]
    bending_axes=['x']

    # Analysis_type= [  'GMNIA','GNA','GNA_Notional_Loads']  
    Analysis_type= [  'GMNIA','GNA','GNA_Notional_Loads']  
    Material_type="50_ksi"
    csv_path = "Column_Results/interaction_results.csv"

    Frame_number=[]
    for column_section_name in column_section_names:
        for story_height in story_heights:
            for no_of_story in No_of_stories:
                for bending_axis in bending_axes:

                    frame_name=check_and_create_new_entries_in_column_config_file(column_section_name,story_height,no_of_story,bending_axis)
                    Frame_number.append(frame_name)


    # Radial_errors=['Error(GNIA<GNA_Notional_Loads)',
    #             'Error(GNIA<GNA)',
    #             'Error(GNA_Notional_Loads<GNA)',
    #             'Error(GNA<GMNIA)',
    #             'Error(GNA_Notional_Loads<GMNIA)',
    #             'Error(GNIA<GMNIA)',
    #             'Error(GMNA<GMNIA)']
    Radial_errors=[
           
            'Error(GNA_Notional_Loads<GNA)',
            'Error(GNA<GMNIA)',
            'Error(GNA_Notional_Loads<GMNIA)',
            
            ]

    if new_analysis_run:
        # Interaction_Plots(
        #     Frame_number=Frame_number,
        #     Analysis_type=Analysis_type,
        #     Material_type=Material_type,
        #     proportional=False,
        #     plot=True)
            
        theta_list = np.linspace(0, 90,91)  
        
        df = write_interaction_results(
            "Column_Results/interaction_results.csv",
            frame_list=Frame_number,
            analysis_list=Analysis_type,
            Material_type=Material_type,
            theta_list=theta_list,
            proportional=False
        )
        # input()
        # plot_theta_vs_del2_over_del1(
        #     csv_path=csv_path,
        #     frame_list=Frame_number,
        #     analyses_to_plot=Analysis_type,
        # )

        # plot_theta_vs_Radial_Errors(
        # csv_path=csv_path,
        # frame_list=Frame_number,
        # errors_to_plot=Radial_errors)

        # parametric_h_min_radial_error(
        #     csv_path=csv_path,
        #     column_section_name="W8X31",
        #     bending_axes="x",
        #     no_of_stories=1,
        #     plot_slenderness=True  
        # )


    else:
        plot_theta_vs_del2_over_del1(
            csv_path=csv_path,
            frame_list=Frame_number,
            analyses_to_plot=Analysis_type,
        )

        plot_theta_vs_Radial_Errors(
        csv_path=csv_path,
        frame_list=Frame_number,
        errors_to_plot=Radial_errors)

        parametric_h_min_radial_error(
            csv_path=csv_path,
            column_section_name="W14X120",
            bending_axes="x",
            no_of_stories=1,
            plot_slenderness=False  
        )

        parametric_slenderness_ratio_vs_min_radial_error(
            csv_path=csv_path)
