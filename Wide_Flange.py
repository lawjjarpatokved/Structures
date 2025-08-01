import libdenavit.section.database.aisc as section
from libdenavit.section.geometric_shape import *
from libdenavit.OpenSees.get_fiber_data import *
from math import pi, ceil
import openseespy.opensees as ops
from Units import *

class WF_Database:
    def __init__(self,Section_name,unit=inch):

        self.section=Section_name
        self.d=section.wide_flange_database[self.section]['d']*unit
        self.tw=section.wide_flange_database[self.section]['tw']*unit
        self.bf=section.wide_flange_database[self.section]['bf']*unit
        self.tf=section.wide_flange_database[self.section]['tf']*unit
        self.A=section.wide_flange_database[self.section]['A']*(unit**2)
        self.Ix=section.wide_flange_database[self.section]['Ix']*(unit**4)
        self.Iy=section.wide_flange_database[self.section]['Iy']*(unit**4)

class I_shape(GeometricShape):
        
    def __init__(self,d,tw,bf,tf,A=None,Ix=None,Iy=None):

        self.d=d
        self.tw=tw
        self.bf=bf
        self.tf=tf
        self._A=A
        self._Ix=Ix
        self._Iy=Iy

        self.num_regions = 10   # Number of regions for the discretization of residual stress
        self.L = float('nan')  # Beam length (span)
        self.TW = float('nan')  # Tributary width
        self.gamma = 0.0

        self.zi = 0.0
        self.zj = 0.0
        self.c = 0.0           # Camber

        self.wD = 0.0           # Dead load (force per unit length)

        self.material_type = 'Steel01'
        self.num_elements = 20
        self.num_fiber = 20
        self.num_steps = 100

        self.nsteps_vol = 30
        self.max_volume = float('nan')
        self.vol_tol = float('nan')
        self.percent_drop = 0.05


    @property
    def A(self):
        if self._A is not None:
            return self._A
        else:
            return 2 * self.bf * self.tf + (self.d - 2 * self.tf) * self.tw

    @A.setter
    def A(self, x):
        self._A = x

    @property
    def Ix(self):
        if self._Ix is not None:
            return self._Ix
        else:
            flange = (self.bf * self.tf**3) / 6
            web = (self.tw * (self.d - 2 * self.tf)**3) / 12
            return 2 * flange + web

    @Ix.setter
    def Ix(self, x):
        self._Ix = x

    @property
    def Iy(self):
        if self._Iy is not None:
            return self._Iy
        else:
            flange = (self.tf * self.bf**3) / 12
            web = ((self.d - 2 * self.tf) * self.tw**3) / 12
            return 2 * flange + web

    @Iy.setter
    def Iy(self, x):
        self._Iy = x

    @property
    def dw(self):
        return self.d-2*self.tf
    
    def Iz(self):
        Iz = (1.0/12)*self.bf*self.d**3 - (1.0/12) * \
            (self.bf-self.tw)*self.dw()**3
        return Iz

    def Sz(self):
        Sz = self.Iz()/(self.d/2)
        return Sz

    def Zz(self):
        Zz = 2*((self.tf*self.bf)*(self.d/2-self.tf/2) +
                (self.dw()/2*self.tw)*(self.dw()/4))
        return Zz

    def Mp(self):
        return self.Fy*self.Zz()

    def C(self):
        return (self.gamma*self.TW*self.L**4)/(pi**4*self.E*self.Iz())
    

    def W_Fiber_section_major_axis(self, sec_tag, material, section_name, Residual_Stress):
        self.sec_tag=sec_tag
        self.Residual_Stress = Residual_Stress

        if self.Residual_Stress:
            self.frc = -0.3 * material.Fy
        else:
            self.frc = 0

        Nfw = ceil(self.dw * (self.num_fiber / self.d))
        Nff = ceil(self.tf * (self.num_fiber / self.d))

        if not self.Residual_Stress:
            ops.section('Fiber', sec_tag)

            ops.uniaxialMaterial('Steel01', material.mat_tag, material.Fy, material.E, material.b)

            # Web patch
            ops.patch('rect', material.mat_tag, Nfw, 1, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)
            # Flange patches
            ops.patch('rect', material.mat_tag, Nff, 1, self.dw / 2, -self.bf / 2, self.d / 2, self.bf / 2)
            ops.patch('rect', material.mat_tag, Nff, 1, -self.d / 2, -self.bf / 2, -self.dw / 2, self.bf / 2)


        else:
            ops.section('Fiber', sec_tag)
            

            frt = -self.frc * (self.bf * self.tf) / (self.bf * self.tf + self.tw * self.dw)

            ops.uniaxialMaterial('Steel01', material.mat_tag + 1, material.Fy, material.E, material.b)
            ops.uniaxialMaterial('InitStressMaterial', material.mat_tag, material.mat_tag + 1, frt)

            ops.patch('rect', material.mat_tag, Nfw, 1, -self.dw / 2, -self.tw / 2, self.dw / 2, self.tw / 2)

            region_width = self.bf / self.num_regions

            for i in range(self.num_regions):
                fri = self.frc + ((i + 0.5) / self.num_regions) * (frt - self.frc)
                matTagi = material.mat_tag + 2 * (i + 1)

                ops.uniaxialMaterial('Steel01', matTagi + 1, material.Fy, material.E, material.b)
                ops.uniaxialMaterial('InitStressMaterial', matTagi, matTagi + 1, fri)

                # Top flange
                ops.patch('rect', matTagi, Nff, 1, self.dw / 2, -region_width / 2, self.d / 2, region_width / 2)
                # Bottom flange
                ops.patch('rect', matTagi, Nff, 1, -self.d / 2, -region_width / 2, -self.dw / 2, region_width / 2)



    def W_Fiber_section_minor_axis(self, sec_tag, material, section_name, Residual_Stress, plot=True):
        self.sec_tag=sec_tag
        self.Residual_Stress = Residual_Stress

        if self.Residual_Stress:
            self.frc = -0.3 * material.Fy
        else:
            self.frc = 0

        Nfw = ceil(self.tw * (self.num_fiber / self.bf))
        Nff = ceil(self.bf * (self.num_fiber / self.bf))

        if not self.Residual_Stress:
            ops.section('Fiber', sec_tag)
            ops.uniaxialMaterial('Steel01', material.mat_tag, material.Fy, material.E, material.b)

            # Web patch
            ops.patch('rect', material.mat_tag, Nfw, 1, -self.tw/2, -self.dw/2, self.tw/2, self.dw/2)

            # Flange patches
            ops.patch('rect', material.mat_tag, Nff, 1, -self.bf/2, -self.d/2, self.bf/2, -self.dw/2)
            ops.patch('rect', material.mat_tag, Nff, 1, -self.bf/2, self.dw/2, self.bf/2, self.d/2)


        else:
            ops.section('Fiber', sec_tag)

            frt = -self.frc * (self.bf * self.tf) / (self.bf * self.tf + self.tw * self.dw)
            ops.uniaxialMaterial('Steel01', material.mat_tag + 1, material.Fy, material.E, material.b)
            ops.uniaxialMaterial('InitStressMaterial', material.mat_tag, material.mat_tag + 1, frt)

            # Web patch
            ops.patch('rect', material.mat_tag, Nfw, 1, -self.tw/2, -self.dw/2, self.tw/2, self.dw/2)

            region_width = self.bf / self.num_regions
            Nff = ceil(region_width * (self.num_fiber / self.bf))

            for i in range(self.num_regions):
                fri = self.frc + ((i + 0.5) / self.num_regions) * (frt - self.frc)
                matTagi = material.mat_tag + 2 * (i + 1)

                ops.uniaxialMaterial('Steel01', matTagi + 1, material.Fy, material.E, material.b)
                ops.uniaxialMaterial('InitStressMaterial', matTagi, matTagi + 1, fri)

                # Top flange segment
                ops.patch('rect', matTagi, Nff, 1, -region_width / 2, self.dw/2, region_width / 2, self.d/2)
                # Bottom flange segment
                ops.patch('rect', matTagi, Nff, 1, -region_width / 2, -self.d/2, region_width / 2, -self.dw/2)



    def plot_fiber_section(sec_tag):
        get_fiber_data(section_tag=sec_tag,plot_fibers=True)


        