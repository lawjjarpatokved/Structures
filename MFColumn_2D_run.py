from Moment_Frame_2D_Main import *
from MFColumn2D import *
from Ziemian_database import Frame_Info,Analysis_Info
from Ziemian_database import convert_dict_items_to_class_attributes
from libdenavit.OpenSees.get_fiber_data import *
from libdenavit.OpenSees import plotting
from Plots import plot_single_bar,line_plot
import os
import seaborn as sns
from typing import List, Optional, Tuple
from Units import load_wind_dirn_data,save_wind_dirn_data,ensure_frame_entry_exists

json_wind_dirn_path="wind_load_dirn_data.json"


def MF_2D_runner(Frame_number,Analysis_type,control_dir='L',lateral_load_scale=1,vertical_load_scale=1,ops_anlaysis='proportional_limit_point'):
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
    results,fail_during_LCA=Frame.run_displacement_controlled_analysis(target_disp=5,plot_defo=False,control_dir=control_dir,
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

def Bar_plot_comparison(Frame_number,Analysis_type):
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

def Interaction_Plots(Frame_number,Analysis_type,proportional=False,plot='False'):
    ##non proportional
    if not proportional:
        code="Non_Proportional"
        for frame_number in Frame_number:
            palette = sns.color_palette("tab10", len(Analysis_type))
            linestyles = ['-', '-', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 1))]

            fig_alr, ax_alr = plt.subplots(figsize=(5, 5))

            # store intersection points for optional later use
            intersections = []   # list of (x, y, color)

            for j, analysis_type in enumerate(Analysis_type):
                ALR_H, ALR_V = [], []

                results, frame_id, fail_during_LCA,Frame = MF_2D_runner(
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
                for i in np.arange(0, 0.2, 0.05):
                    print(f"Running vertical load scale {i:.2f}")
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
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

                # --- Sweep vertical loads (0.8–1.0) ---
                for i in np.arange(0.2, 0.9, 0.1):
                    print(f"Running vertical load scale {i:.2f}")
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
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

                for i in np.arange(0.9, 1.01, 0.01):
                    print(f"Running vertical load scale {i:.2f}")
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
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


    else:
    ### proportional
        code="Proportional"
        for frame_number in Frame_number:
            palette = sns.color_palette("tab10", len(Analysis_type))
            linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1)), (0, (5, 1))]

            fig_alr, ax_alr = plt.subplots(figsize=(5, 5))

            # NEW: store intersections if you want them later (optional)
            intersections = []

            for j, analysis_type in enumerate(Analysis_type):
                ALR_H, ALR_V = [], []

                # --- Base case: vertical-controlled analysis ---
                results, frame_id, fail_during_LCA,Frame = MF_2D_runner(
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
                for i in np.arange(0, 0.1, 0.05):
                    
                    print(f"Running vertical load scale {i:.3f}")
                    
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
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
                    
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point'
                    )

                    if fail_during_LCA or results.maximum_load_ratio_at_limit_point < 0.01:
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
                for i in np.arange(1.0, 6.0, 0.1):
                    
                    print(f'{analysis_type}')
                    print(f"Running vertical load scale {i:.2f}")
                    input()

    
                    
                    results, frame_id, fail_during_LCA,_ = MF_2D_runner(
                        Frame_number=frame_number,
                        Analysis_type=analysis_type,
                        vertical_load_scale=i * ALR_V_max,
                        control_dir='L',
                        ops_anlaysis='proportional_limit_point'
                    )

                    if fail_during_LCA or results.maximum_load_ratio_at_limit_point < 0.01:
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
        ax_alr.legend()
        fig_alr.tight_layout()
        fig_alr.savefig(
            os.path.join(frame_id, f'{code}_ALR_H_vs_ALR_V_Frame{frame_number}.png'),
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
        if hcol not in df.columns:
            df[hcol] = np.nan
        if vcol not in df.columns:
            df[vcol] = np.nan

    # --- Index for updates ---
    df.set_index(["Frame", "Theta"], inplace=True)

    # --- Compute + fill only the analyses requested ---
    for frame in frame_list:
        frame = str(frame)
        print(f"\n=== Frame: {frame} ===")

        curves = {}
        for analysis in analysis_list:
            print(f"Running {analysis} for {frame}")
            ALR_H, ALR_V,duplicate_frame = Interaction_Plots(
                Frame_number=[frame],
                Analysis_type=[analysis],
                proportional=proportional,
                plot=False
            )
            curves[analysis] = (ALR_H, ALR_V)
        print(theta_list)
        for analysis in analysis_list:
            for theta in theta_list:
                ALR_H, ALR_V = curves[analysis]

                pt = intersection_with_ray_from_origin(ALR_H, ALR_V, theta)
                print(theta)
                print(pt)
                # input()
                vertical_load_scale = pt[1] if pt is not None else None
                if vertical_load_scale==0:
                    ult_lat_load_for_V0=pt[0]
                    print(vertical_load_scale)
                    print(ult_lat_load_for_V0)
                    input()
                lateral_load_scale = ult_lat_load_for_V0*0.0001 if pt[0]==0 else pt[0] if pt is not None else None
                # lateral_load_scale = pt[0]
                del2_over_del1=duplicate_frame.get_del2_over_del1(vertical_load_scale=vertical_load_scale, lateral_load_scale=lateral_load_scale)


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
    os.makedirs(os.path.dirname(csv_path) or ".", exist_ok=True)
    df.to_csv(csv_path, index=False)

    return df


def return_radial_error_betn_analyses(df: pd.DataFrame, analysis1: str, analysis2: str) -> pd.Series:
    """
    Calculate the radial error between two analyses for each row in the DataFrame.
    The radial error is defined as the Manhattan distance between the two analyses
    """
    ana1_ALR_H = f"{analysis1}_ALR_H"
    ana1_ALR_V = f"{analysis1}_ALR_V"
    ana2_ALR_H = f"{analysis2}_ALR_H"
    ana2_ALR_V = f"{analysis2}_ALR_V"

    if ana1_ALR_H not in df.columns or ana2_ALR_H not in df.columns or ana1_ALR_V not in df.columns or ana2_ALR_V not in df.columns:
        raise ValueError(f"Columns '{ana1_ALR_H}' and/or '{ana2_ALR_H}' and/or '{ana1_ALR_V}' and/or '{ana2_ALR_V}' not found in DataFrame.")

    val1_H = pd.to_numeric(df[ana1_ALR_H], errors="coerce")
    val2_H = pd.to_numeric(df[ana2_ALR_H], errors="coerce")
    val1_V = pd.to_numeric(df[ana1_ALR_V], errors="coerce")
    val2_V = pd.to_numeric(df[ana2_ALR_V], errors="coerce")

    radial_error_H = (val2_H - val1_H) 
    radial_error_V = (val2_V - val1_V)
    radial_error = (radial_error_H + radial_error_V)
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
    else:
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
    else:
        plt.close()

if __name__ == "__main__":
    new_analysis_run=True
    Frame_number= ['Trial_Col_2']      # 'SP36H'  ,  'UP36H'  ,  'SP36L'  ,  'UP36L'
    Analysis_type= ['GMNA'  ,  'GMNIA'  ,'GNA', 'GNIA', 'GNA_Notional_Loads']

    Analysis_type= [  'GMNIA','GNA','GNA_Notional_Loads']
    csv_path = "Column_Results/interaction_results.csv"


    Radial_errors=['Error(GNIA<GNA_Notional_Loads)',
                'Error(GNIA<GNA)',
                'Error(GNA_Notional_Loads<GNA)',
                'Error(GNA<GMNIA)',
                'Error(GNA_Notional_Loads<GMNIA)',
                'Error(GNIA<GMNIA)',
                'Error(GMNA<GMNIA)']

    if new_analysis_run:
        # Interaction_Plots(
        #     Frame_number=Frame_number,
        #     Analysis_type=Analysis_type,
        #     proportional=False,
        #     plot=True)
        # input('Plots Created')    
        theta_list = np.linspace(0, 90,91)  
        
        df = write_interaction_results(
            "Column_Results/interaction_results.csv",
            frame_list=Frame_number,
            analysis_list=Analysis_type,
            theta_list=theta_list,
            proportional=False
        )
        
        plot_theta_vs_del2_over_del1(
            csv_path=csv_path,
            frame_list=Frame_number,
            analyses_to_plot=Analysis_type,
        )

        # plot_theta_vs_Radial_Errors(
        # csv_path=csv_path,
        # frame_list=Frame_number,
        # errors_to_plot=Radial_errors)


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