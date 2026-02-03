from Moment_Frame_2D_Main import *
from Ziemian_database import Frame_Info,Analysis_Info
from Ziemian_database import convert_dict_items_to_class_attributes
from libdenavit.OpenSees.get_fiber_data import *
from libdenavit.OpenSees import plotting
import gc
from Plots import plot_single_bar,line_plot
import os
import seaborn as sns
import opsvis

def MF_2D_runner(Frame_number,Analysis_type,control_dir='L',lateral_load_scale=1,vertical_load_scale=1,ops_anlaysis='proportional_limit_point'):
    # try:
    Frame_dict=Frame_Info[str(Frame_number)]
    Frame_details=convert_dict_items_to_class_attributes(Frame_dict)
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
                    geometric_imperfection_ratio=Frame_details.geometric_imperfection_ratio,
                    nip=3,
                    mat_type='Steel01',
                    wind_load_dirn=wind_load_dirn,
                    Leaning_column=Frame_details.Leaning_column,
                    Leaning_column_offset=Frame_details.Leaning_column_offset,
                    Leaning_column_floor_load=Frame_details.Leaning_column_floor_load,
                    Leaning_column_roof_load=Frame_details.Leaning_column_roof_load)

    Frame.generate_Nodes_and_Element_Connectivity()
    Frame.create_distorted_nodes_and_element_connectivity()
    Frame.build_ops_model()
    results,fail_during_LCA=Frame.run_displacement_controlled_analysis(target_disp=5,plot_defo=False,control_dir=control_dir,
                lateral_load_scale=lateral_load_scale,vertical_load_scale=vertical_load_scale,analysis=ops_anlaysis)
    
    # Frame.plot_model()
    return results,Frame_details.Frame_id,fail_during_LCA

    # finally:
    #     # --- CLEANUP ---
    #     ops.wipe()          # clear OpenSees model
    #     if 'Frame' in locals():
    #         Frame.__dict__.clear()
    #         del Frame          # delete Python object
    #     gc.collect()

def Bar_plot_comparison(Frame_number,Analysis_type):
    for frame_number in Frame_number:
        max_load_ratio = []
        analysis_type_labels = []

        for analysis_type in Analysis_type:
            print(f'Running {analysis_type}')
            analysis_type_labels.append(analysis_type.replace("_", " "))

            # --- Run analysis ---
            results, frame_id, fail_during_LCA = MF_2D_runner(
                Frame_number=frame_number,
                Analysis_type=analysis_type,
                control_dir='L',
                ops_anlaysis='proportional_limit_point'
            )
            # --- Create folder structure ---
            # frame_id is usually like "Frame_1"
            analysis_folder = os.path.join(frame_id, analysis_type)
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
        os.makedirs(frame_id, exist_ok=True)
        barplot_filename = os.path.join(frame_id, f"{frame_number}_Barplot.png")
        plot_single_bar(max_load_ratio, filename=barplot_filename,
                        x_labels=analysis_type_labels,
                        title=f'{frame_number}: Comparison of Load Ratio',
                        ylabel='Load Ratio',title_fontsize=16, label_fontsize=14, tick_fontsize=14, legend_fontsize=12,value_fontsize=12)

def Interaction_Plots(Frame_number,Analysis_type,proportional=False):
    ##non proportional
    if not proportional:
        for frame_number in Frame_number:
            palette = sns.color_palette("tab10", len(Analysis_type))
            linestyles = ['-', '-', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 1))]

            fig_alr, ax_alr = plt.subplots(figsize=(5, 5))

            # store intersection points for optional later use
            intersections = []   # list of (x, y, color)

            for j, analysis_type in enumerate(Analysis_type):
                ALR_H, ALR_V = [], []

                results, frame_id, fail_during_LCA = MF_2D_runner(
                    Frame_number=frame_number,
                    Analysis_type=analysis_type,
                    lateral_load_scale=0,
                    control_dir='V',
                    ops_anlaysis='proportional_limit_point'
                )

                analysis_folder = os.path.join(frame_id, analysis_type)
                os.makedirs(analysis_folder, exist_ok=True)

                # --- Initialize ALR values ---
                ALR_V_max = results.maximum_load_ratio_at_limit_point
                ALR_V.append(ALR_V_max)
                ALR_H.append(0)

                # --- Plot PMM interaction for base case ---
                fig_pmm, ax_pmm = plotting.plot_PMM_Interaction_values(
                    results.P_M_M_interaction_all_elements[-1],
                    show=False
                )
                fig_pmm.savefig(
                    os.path.join(
                        analysis_folder,
                        f"PMM_{analysis_type}_ALRH_{ALR_H[0]:.2f}_ALRV_{ALR_V[0]:.2f}.png"
                    ),
                    dpi=600
                )
                plt.close(fig_pmm)

                # --- Sweep vertical loads (0–0.8) ---
                for i in np.arange(0, 0.5, 0.05):
                    print(f"Running vertical load scale {i:.2f}")
                    results, frame_id, fail_during_LCA = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='non_proportional_limit_point'
                    )
                    if fail_during_LCA:
                        break
                    else:
                        ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                        ALR_V.insert(-1, i * ALR_V_max)

                        if int(i * 10) % 2 == 0:
                            fig_pmm_i, ax_pmm_i = plotting.plot_PMM_Interaction_values(
                                results.P_M_M_interaction_all_elements[-1],
                                show=False
                            )
                            fig_pmm_i.savefig(
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                ),
                                dpi=600
                            )
                            plt.close(fig_pmm_i)

                # # --- Sweep vertical loads (0.8–1.0) ---
                # for i in np.arange(0.2, 0.9, 0.1):
                #     print(f"Running vertical load scale {i:.2f}")
                #     results, frame_id, fail_during_LCA = MF_2D_runner(
                #         Frame_number=frame_number,
                #         Analysis_type=analysis_type,
                #         vertical_load_scale=i * ALR_V_max,
                #         control_dir='L',
                #         ops_anlaysis='non_proportional_limit_point'
                #     )
                #     if fail_during_LCA:
                #         break
                #     else:
                #         ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                #         ALR_V.insert(-1, i * ALR_V_max)

                #         if int(i * 10) % 2 == 0:
                #             fig_pmm_i, ax_pmm_i = plotting.plot_PMM_Interaction_values(
                #                 results.P_M_M_interaction_all_elements[-1],
                #                 show=False
                #             )
                #             fig_pmm_i.savefig(
                #                 os.path.join(
                #                     analysis_folder,
                #                     f"PMM_{analysis_type}_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                #                 ),
                #                 dpi=600
                #             )
                #             plt.close(fig_pmm_i)

                # for i in np.arange(0.9, 1.01, 0.01):
                #     print(f"Running vertical load scale {i:.2f}")
                #     results, frame_id, fail_during_LCA = MF_2D_runner(
                #         Frame_number=frame_number,
                #         Analysis_type=analysis_type,
                #         vertical_load_scale=i * ALR_V_max,
                #         control_dir='L',
                #         ops_anlaysis='non_proportional_limit_point'
                #     )
                #     if fail_during_LCA:
                #         break
                #     else:
                #         ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                #         ALR_V.insert(-1, i * ALR_V_max)

                #         if int(i * 10) % 2 == 0:
                #             fig_pmm_i, ax_pmm_i = plotting.plot_PMM_Interaction_values(
                #                 results.P_M_M_interaction_all_elements[-1],
                #                 show=False
                #             )
                #             fig_pmm_i.savefig(
                #                 os.path.join(
                #                     analysis_folder,
                #                     f"PMM_{analysis_type}_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                #                 ),
                #                 dpi=600
                #             )
                #             plt.close(fig_pmm_i)

                # Convert to arrays for intersection math
                ALR_H_arr = np.array(ALR_H, dtype=float)
                ALR_V_arr = np.array(ALR_V, dtype=float)

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

                # ---------- NEW: find intersection with line y = x ----------
                diff = ALR_V_arr - ALR_H_arr
                # indices where sign changes: diff[i] * diff[i+1] <= 0
                sign_change_idx = np.where(diff[:-1] * diff[1:] <= 0)[0]

                if sign_change_idx.size > 0:
                    k = sign_change_idx[0]
                    x1, x2 = ALR_H_arr[k], ALR_H_arr[k+1]
                    y1, y2 = ALR_V_arr[k], ALR_V_arr[k+1]

                    # parametric interpolation to solve for x where y=x
                    # segment: (x(t), y(t)) = (x1 + t*(x2-x1), y1 + t*(y2-y1))
                    # find t such that x(t) = y(t):
                    denom = (x2 - x1) - (y2 - y1)
                    if abs(denom) > 1e-9:
                        t = (y1 - x1) / denom
                        t = np.clip(t, 0.0, 1.0)
                        x_cross = x1 + t * (x2 - x1)
                        y_cross = y1 + t * (y2 - y1)
                    else:
                        # nearly parallel; just use midpoint
                        x_cross = 0.5 * (x1 + x2)
                        y_cross = 0.5 * (y1 + y2)

                    intersections.append((x_cross, y_cross, palette[j % len(palette)]))
                    # plot small highlight marker immediately
                    ax_alr.scatter(
                        x_cross, y_cross,
                        s=15,
                        facecolors='none',
                        edgecolors=palette[j % len(palette)],
                        linewidths=0.8,
                        zorder=5
                    )

            # --- After all curves: add diagonal reference line and finalize ---
            # 45° line from (0,0) to (1,1)
            ax_alr.plot([0, 1], [0, 1],color='black',linestyle='--',linewidth=0.8,alpha=0.7,label='ALR_H = ALR_V')

            ax_alr.set_xlim(left=0)
            ax_alr.set_ylim(bottom=0)

            ax_alr.set_title(f'Frame {frame_number}: ALR_V vs ALR_H')
            ax_alr.legend()
            fig_alr.tight_layout()
            fig_alr.savefig(os.path.join(frame_id, f'Non_prop_ALR_H_vs_ALR_V_Frame{frame_number}.png'), dpi=600,transparent=False)
            plt.close(fig_alr)

            # plotting.plot_sfd()
            # plotting.plot_bmd()
            # plotting.plot_afd(scale=0.001)




    else:
    ### proportional

        for frame_number in Frame_number:
            palette = sns.color_palette("tab10", len(Analysis_type))
            linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 1))]

            fig_alr, ax_alr = plt.subplots(figsize=(5, 5))

            # NEW: store intersections if you want them later (optional)
            intersections = []

            for j, analysis_type in enumerate(Analysis_type):
                ALR_H, ALR_V = [], []

                # --- Base case: vertical-controlled analysis ---
                results, frame_id, fail_during_LCA = MF_2D_runner(
                    Frame_number=frame_number,
                    Analysis_type=analysis_type,
                    lateral_load_scale=0,
                    control_dir='V',
                    ops_anlaysis='proportional_limit_point'
                )

                analysis_folder = os.path.join(frame_id, analysis_type)
                os.makedirs(analysis_folder, exist_ok=True)

                # --- Initialize ALR values ---
                ALR_V_max = results.maximum_load_ratio_at_limit_point
                ALR_V.append(ALR_V_max)
                ALR_H.append(0)

                # --- Plot PMM interaction for base case ---
                fig_pmm, ax_pmm = plotting.plot_PMM_Interaction_values(
                    results.P_M_M_interaction_all_elements[-1],
                    show=False
                )
                fig_pmm.savefig(
                    os.path.join(
                        analysis_folder,
                        f"PMM_{analysis_type}_ALRH_{ALR_H[0]:.2f}_ALRV_{ALR_V[0]:.2f}.png"
                    ),
                    dpi=600
                )
                plt.close(fig_pmm)

                # --- Sweep vertical loads: fine steps near 0 ---
                for i in np.arange(0, 0.1, 0.005):
                    print(f"Running vertical load scale {i:.3f}")
                    results, frame_id, fail_during_LCA = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
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

                        # --- Plot PMM every 2nd step ---
                        if int(i * 10) % 2 == 0:
                            fig_pmm_i, ax_pmm_i = plotting.plot_PMM_Interaction_values(
                                results.P_M_M_interaction_all_elements[-1],
                                show=False
                            )
                            fig_pmm_i.savefig(
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_proportional_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                ),
                                dpi=600
                            )
                            plt.close(fig_pmm_i)

                # --- Sweep vertical loads: 0.1 to 1.0 ---
                for i in np.arange(0.1, 1.0, 0.1):
                    print(f"Running vertical load scale {i:.2f}")
                    results, frame_id, fail_during_LCA = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point'
                    )

                    if fail_during_LCA and results.maximum_load_ratio_at_limit_point < 0.01:
                        break
                    else:
                        ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                        ALR_V.insert(-1, i * ALR_V_max * results.maximum_load_ratio_at_limit_point)

                        if int(i * 10) % 2 == 0:
                            fig_pmm_i, ax_pmm_i = plotting.plot_PMM_Interaction_values(
                                results.P_M_M_interaction_all_elements[-1],
                                show=False
                            )
                            fig_pmm_i.savefig(
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_proportional_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                ),
                                dpi=600
                            )
                            plt.close(fig_pmm_i)

                # --- Sweep vertical loads: 1.0 to 4.0 ---
                for i in np.arange(1.0, 4.0, 0.02):
                    print(f"Running vertical load scale {i:.2f}")
                    results, frame_id, fail_during_LCA = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point'
                    )

                    if fail_during_LCA and results.maximum_load_ratio_at_limit_point < 0.01:
                        break
                    else:
                        ALR_H.insert(-1, results.maximum_load_ratio_at_limit_point)
                        ALR_V.insert(-1, i * ALR_V_max * results.maximum_load_ratio_at_limit_point)

                        if int(i * 10) % 2 == 0:
                            fig_pmm_i, ax_pmm_i = plotting.plot_PMM_Interaction_values(
                                results.P_M_M_interaction_all_elements[-1],
                                show=False
                            )
                            fig_pmm_i.savefig(
                                os.path.join(
                                    analysis_folder,
                                    f"PMM_{analysis_type}_proportional_ALRH_{ALR_H[-2]:.2f}_ALRV_{ALR_V[-2]:.2f}.png"
                                ),
                                dpi=600
                            )
                            plt.close(fig_pmm_i)

                # ============================================================
                #   PLOTTING ALR_H vs ALR_V FOR THIS ANALYSIS TYPE
                # ============================================================

                ALR_H_arr = np.array(ALR_H, dtype=float)
                ALR_V_arr = np.array(ALR_V, dtype=float)

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

                # ---------- Intersection with line y = x ----------
                diff = ALR_V_arr - ALR_H_arr
                sign_change_idx = np.where(diff[:-1] * diff[1:] <= 0)[0]

                if sign_change_idx.size > 0:
                    k = sign_change_idx[0]
                    x1, x2 = ALR_H_arr[k], ALR_H_arr[k + 1]
                    y1, y2 = ALR_V_arr[k], ALR_V_arr[k + 1]

                    denom = (x2 - x1) - (y2 - y1)
                    if abs(denom) > 1e-9:
                        t = (y1 - x1) / denom
                        t = np.clip(t, 0.0, 1.0)
                        x_cross = x1 + t * (x2 - x1)
                        y_cross = y1 + t * (y2 - y1)
                    else:
                        x_cross = 0.5 * (x1 + x2)
                        y_cross = 0.5 * (y1 + y2)

                    intersections.append((x_cross, y_cross, palette[j % len(palette)]))

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

            ax_alr.set_title(f'Frame {frame_number}: ALR_V vs ALR_H (Proportional)')
            ax_alr.legend()
            fig_alr.tight_layout()
            fig_alr.savefig(
                os.path.join(frame_id, f'Proportional_ALR_H_vs_ALR_V_Frame{frame_number}.png'),
                dpi=600
            )
            plt.close(fig_alr)

            # plotting.plot_sfd()
            # plotting.plot_bmd()
            # plotting.plot_afd(scale=0.001)



Frame_number= ['SF36H']      # 'SP36H'  ,  'UP36H'  ,  'SP36L'  ,  'UP36L'
# Analysis_type= ['GMNA'  ,  'GMNIA'  ,'GNA', 'GNIA', 'GNA_Notional_Loads']
Analysis_type= [ 'GMNIA']


# Bar_plot_comparison(Frame_number=Frame_number,Analysis_type=Analysis_type)
Interaction_Plots(Frame_number=Frame_number,Analysis_type=Analysis_type,proportional=False)

