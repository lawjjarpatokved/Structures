import numpy as np
import math as math
import pandas as pd
import openseespy.opensees as ops
import os
import opsvis as opsv
from libdenavit.OpenSees.plotting import *
from libdenavit.OpenSees.get_fiber_data import *
from libdenavit.cross_section_2d import *
from Units import *
from math import pi, ceil
from libdenavit.section.wide_flange import *
from Moment_Frame_2D_Main import Steel_Material


#################################################
density_of_steel=7850*kg/(m**3)
g=9.81*m/(sec**2)
E=29000*ksi
G=77221*Mpa
fy=36*ksi
Hk = 0.001*E           # Kinematic hardening modulus

beam_section_tag=1
Steel=Steel_Material(1,E=E,fy=fy,G=G,Hk=Hk,density=density_of_steel)
beam_section_name="W27X84"
beam_data = wf_Database(beam_section_name)

Mp=Steel.fy*beam_data.Zx
Initial_slope=Steel.E*beam_data.Ix


def Moment_curvature_analysis(P_value, residual_stress, axis):
    frc = -0.3 * Steel.fy if residual_stress else 0
    beam = I_shape(beam_data.d, beam_data.tw, beam_data.bf, beam_data.tf,
                   fy=Steel.fy, E=Steel.E, Hk=Steel.Hk,
                   A=beam_data.A, Ix=beam_data.Ix, Iy=beam_data.Iy)
    member = CrossSection2d(beam, axis=axis)
    results = member.run_ops_analysis(
        analysis_type='nonproportional_limit_point',
        **{
            'section_id': beam_section_tag,
            'section_args': [Steel.mat_tag, 'Steel01', 20, 20, frc],
            'section_kwargs': {},
            'P': P_value
        }
    )
    return results


def axial_load_comparison(P_values=(0, 500, 1000, 1500, 2000), residual_stress=True, axis='x'):
    plt.figure()
    for P in P_values:
        results = Moment_curvature_analysis(P_value=P, residual_stress=residual_stress, axis=axis)
        label = f"P={results.applied_axial_load[-1]:.0f}"
        plt.plot(results.curvatureX, results.maximum_abs_moment, label=label)

    # plt.axhline(y=Mp, linestyle='--', linewidth=1.2, label=f"M_p = {Mp:.3g}")
    # kappa_end = Mp / Initial_slope
    # plt.plot([0, kappa_end], [0, Mp], linestyle=':', linewidth=1.2, label=f"slope = {Initial_slope:.3g}")

    rs_text = "with RS" if residual_stress else "without RS"
    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title(f"Moment–Curvature vs Axial Load ({rs_text})")
    plt.grid(True)
    plt.legend()
    plt.show()

def comparison_with_EI_and_FyZ(P=0, residual_stress=False, axis='x'):
    plt.figure()
    results = Moment_curvature_analysis(P_value=P, residual_stress=residual_stress, axis=axis)
    label = f"P={results.applied_axial_load[-1]:.0f}"
    plt.plot(results.curvatureX, results.maximum_abs_moment, label=label)

    plt.axhline(y=Mp, linestyle='--', linewidth=1.2, label=f"M_p = {Mp:.3g}")
    kappa_end = Mp / Initial_slope
    plt.plot([0, kappa_end], [0, Mp], linestyle=':', linewidth=1.2, label=f"slope = {Initial_slope:.3g}")

    rs_text = "with RS" if residual_stress else "without RS"
    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title(f"Moment–Curvature vs Axial Load ({rs_text})")
    plt.grid(True)
    plt.legend()
    plt.show()

def residual_stress_comparison():
    P_values = [0,  1000]
    rs_flags = [False, True]  # without RS, with RS
    styles = {False: {'linestyle': '-',  'linewidth': 1.4},
            True:  {'linestyle': '--', 'linewidth': 1.4}}
    plt.figure()

    for P in P_values:
        for rs in rs_flags:
            results = Moment_curvature_analysis(P_value=P, residual_stress=rs, axis='x')
            label = f"P={results.applied_axial_load[-1]:.0f}, RS={'Yes' if rs else 'No'}"
            plt.plot(results.curvatureX, results.maximum_abs_moment, label=label, **styles[rs])

    # reference lines
    # plt.axhline(y=Mp, linestyle='--', linewidth=1.2, label=f"M_p = {Mp:.3g}")
    # kappa_end = Mp / Initial_slope
    # plt.plot([0, kappa_end], [0, Mp], linestyle=':', linewidth=1.2, label=f"slope = {Initial_slope:.3g}")

    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title("Moment–Curvature Curves (with/without Residual Stress)")
    plt.grid(True)
    plt.legend()
    plt.show()

def comparison_of_major_and_minor_axes(P=0, residual_stress=False, axes=('x','y')):
    plt.figure()
    for ax in axes:
        res = Moment_curvature_analysis(P_value=P, residual_stress=residual_stress, axis=ax)
        curv = res.curvatureX if ax == 'x' else res.curvatureY
        plt.plot(curv, res.maximum_abs_moment, label=f"P={res.applied_axial_load[-1]:.0f}, axis={ax}")

        Mp_axis = Steel.fy * (beam_data.Zx if ax == 'x' else beam_data.Zy)
        EI_axis = Steel.E * (beam_data.Ix if ax == 'x' else beam_data.Iy)
        plt.axhline(y=Mp_axis, linestyle='--', linewidth=1.2, label=f"M_p,{ax} = {Mp_axis:.3g}")
        kappa_end = Mp_axis / EI_axis
        plt.plot([0, kappa_end], [0, Mp_axis], linestyle=':', linewidth=1.2, label=f"slope EI_{ax} = {EI_axis:.3g}")

    rs_text = "with RS" if residual_stress else "without RS"
    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title(f"Moment–Curvature (axes={list(axes)}, {rs_text}, P={P})")
    plt.grid(True)
    plt.legend()
    plt.show()

def comparison_of_2d_and_3d_section(P=0, residual_stress=False, axes=('x','y',None)):
    plt.figure()
    for ax in axes:
        res = Moment_curvature_analysis(P_value=P, residual_stress=residual_stress, axis=ax)
        curv = res.curvatureX if ax == 'x' or ax== None else res.curvatureY
        plt.plot(curv, res.maximum_abs_moment, label=f"P={res.applied_axial_load[-1]:.0f}, axis={ax}")

    rs_text = "with RS" if residual_stress else "without RS"
    plt.xlabel("Curvature κ")
    plt.ylabel("Moment M")
    plt.title(f"Moment–Curvature (axes={list(axes)}, {rs_text}, P={P})")
    plt.grid(True)
    plt.legend()
    plt.show()

residual_stress_comparison()
axial_load_comparison()
comparison_with_EI_and_FyZ()
comparison_of_major_and_minor_axes()
comparison_of_2d_and_3d_section()