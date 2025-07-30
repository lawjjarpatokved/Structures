import libdenavit.section.database.aisc as section
from libdenavit.section.geometric_shape import *
from libdenavit.OpenSees.get_fiber_data import *
from math import pi, ceil
import openseespy.opensees as ops
import matplotlib.pyplot as plt
import opsvis as opsv
from Units import *

# class WF_Database:
#     def __init__(self,Section_name):

#         self.section=Section_name
#         self.d=section.wide_flange_database[self.section]['d']*inch
#         self.tw=section.wide_flange_database[self.section]['tw']*inch
#         self.bf=section.wide_flange_database[self.section]['bf']*inch
#         self.tf=section.wide_flange_database[self.section]['tf']*inch
#         self.A=section.wide_flange_database[self.section]['A']*(inch**2)
#         self.Ix=section.wide_flange_database[self.section]['Ix']*(inch**4)
#         self.Iy=section.wide_flange_database[self.section]['Iy']*(inch**4)

class I_shape(GeometricShape):
        
    def __init__(self,Section_name,Fy,E,Hk):

        self.section=Section_name
        self.d=section.wide_flange_database[self.section]['d']*inch
        self.tw=section.wide_flange_database[self.section]['tw']*inch
        self.bf=section.wide_flange_database[self.section]['bf']*inch
        self.tf=section.wide_flange_database[self.section]['tf']*inch
        self.A=section.wide_flange_database[self.section]['A']*(inch**2)
        self.Ix=section.wide_flange_database[self.section]['Ix']*(inch**4)
        self.Iy=section.wide_flange_database[self.section]['Iy']*(inch**4)
        self.E=E
        self.Fy=Fy
        self.Hk = Hk            # Kinematic hardening modulus
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

        self.dw = self.d - 2*self.tf



    def dw(self):
        dw = self.d-2*self.tf
        return dw

    def A(self):
        A = 2*self.bf*self.tf + (self.d-2*self.tf)*self.tw
        return A

    def Iz(self):
        Iz = (1.0/12)*self.bf*self.d**3 - (1.0/12) * \
            (self.bf-self.tw)*self.dw()**3
        return Iz

    def Sz(self):
        Sz = self.Iz()/(self.d/2)
        return Sz

    def Zz(self):
        Zz = 2*((self.tf*self.bf)*(self.d/2-self.tf/2) +
                (self.dw()/2*self.tw)*(self.dw()/4));
        return Zz

    def Mp(self):
        return self.Fy*self.Zz()

    def C(self):
        return (self.gamma*self.TW*self.L**4)/(pi**4*self.E*self.Iz())
    

    def W_Fiber_section_major_axis(self, sec_tag, matTag,section_name,Residual_Stress,plot=True,):
        self.Residual_Stress=Residual_Stress
                # Maximum compressive residual stress (Lehigh pattern)
        if self.Residual_Stress == True:
            self.frc=-0.3*self.Fy
        else:
            self.frc=0
        """
        Creates W-Section based on nominal dimension and generates
        fibers over it   
        """
        Nfw = ceil(self.dw*(self.num_fiber/self.d))
        Nff = ceil(self.tf*(self.num_fiber/self.d))


        if not self.Residual_Stress:

            ops.section('Fiber', sec_tag)
            self.b = self.Hk/(self.E+self.Hk)
            ops.uniaxialMaterial('Steel01', matTag, self.Fy, self.E, self.b)
            ## web
            ops.patch('rect', matTag, Nfw,1,-self.dw/2, -self.tw/2, self.dw/2, self.tw/2)
            ## flanges
            ops.patch('rect', matTag, Nff,1,self.dw/2, -self.bf/2,self.d/2, self.bf/2)
            
            ops.patch('rect', matTag, Nff,1,-self.d/2, -self.bf/2, -self.dw/2, self.bf/2)
            
            
            if plot:

                fib_sec = [

                        ['section', 'Fiber', sec_tag, '-GJ', 1.0e4],

                        ["patch","rect", matTag, Nfw,1,-self.dw/2, -self.tw/2, self.dw/2, self.tw/2],
                        ["patch","rect",matTag, Nff,1,self.dw/2, -self.bf/2,self.d/2, self.bf/2],
                        ["patch","rect", matTag, Nff,1,-self.d/2, -self.bf/2, -self.dw/2, self.bf/2]

                        ]
                
                matcolor = [ 'gold']*10000
                opsv.plot_fiber_section(fib_sec, matcolor=matcolor)
                plt.axis('equal')
                plt.title(section_name)
                plt.show()

        else:
            ops.section('Fiber', sec_tag)


            frt = -self.frc*(self.bf*self.tf) / \
                (self.bf*self.tf+self.tw*self.dw)
            # Define web fibers
            self.b = self.Hk/(self.E+self.Hk)
            ops.uniaxialMaterial('Steel01', matTag+1, self.Fy, self.E, self.b)
            ops.uniaxialMaterial('InitStressMaterial',
                                    matTag, matTag+1, frt)
            ops.patch('rect', matTag, Nfw, 1, -self.dw/2, -self.tw/2, self.dw/2, self.tw/2)

            # Define flange fibers
            region_width = self.bf/self.num_regions
            for i in range(self.num_regions):
                fri = self.frc + ((i+0.5)/self.num_regions)*(frt-self.frc)

                matTagi = matTag+2*(i+1)

                ops.uniaxialMaterial(
                    'Steel01', matTagi+1, self.Fy, self.E, self.b)
                ops.uniaxialMaterial(
                    'InitStressMaterial', matTagi, matTagi+1, fri)

                ops.patch('rect', matTagi, Nff, 1, self.dw/2, -
                          region_width/2, self.d/2, region_width/2)
                ops.patch('rect', matTagi, Nff, 1, -self.d/2, -
                          region_width/2, -self.dw/2, region_width/2)
                
            if plot:

                fib_sec = [

                        ['section', 'Fiber', sec_tag, '-GJ', 1.0e4],

                        ["patch","rect", matTag,Nfw, 1, -self.dw/2, -self.tw/2, self.dw/2, self.tw/2],
                        ["patch","rect",matTagi, Nff, 1, self.dw/2, -region_width/2, self.d/2, region_width/2],
                        ["patch","rect", matTagi, Nff, 1, -self.d/2, -region_width/2, -self.dw/2, region_width/2]

                        ]
                
                matcolor = [ 'gold']*10000
                opsv.plot_fiber_section(fib_sec, matcolor=matcolor)
                plt.axis('equal')
                plt.title(section_name)
                plt.show()


    def W_Fiber_section_minor_axis(self, sec_tag, matTag,section_name,Residual_Stress,plot=True,):
        self.Residual_Stress=Residual_Stress
                # Maximum compressive residual stress (Lehigh pattern)
        if self.Residual_Stress == True:
            self.frc=-0.3*self.Fy
        else:
            self.frc=0
        """
        Creates W-Section based on nominal dimension and generates
        fibers over it   
        """
        Nfw = ceil(self.tw*(self.num_fiber/self.bf))
        Nff = ceil(self.bf*(self.num_fiber/self.bf))
        


        if not self.Residual_Stress:

            ops.section('Fiber', sec_tag)
            self.b = self.Hk/(self.E+self.Hk)
            ops.uniaxialMaterial('Steel01', matTag, self.Fy, self.E, self.b)
            ## web
            ops.patch('rect', matTag, Nfw,1,-self.tw/2,-self.dw/2, self.tw/2,self.dw/2)
            ## flanges
            ops.patch('rect', matTag, Nff, 1, -self.bf/2,-self.d/2, self.bf/2,-self.dw/2)
            
            ops.patch('rect', matTag, Nff, 1,-self.bf/2,self.dw/2, self.bf/2,self.d/2)
            
            
            if plot:

                fib_sec = [

                        ['section', 'Fiber', sec_tag, '-GJ', 1.0e4],

                        ["patch","rect", matTag,Nfw,1,-self.tw/2,-self.dw/2, self.tw/2,self.dw/2],
                        ["patch","rect",matTag, Nff, 1, -self.bf/2,-self.d/2, self.bf/2,-self.dw/2],
                        ["patch","rect", matTag, Nff, 1,-self.bf/2,self.dw/2, self.bf/2,self.d/2]
                ]
                
                matcolor = [ 'gold']*10000
                opsv.plot_fiber_section(fib_sec, matcolor=matcolor)
                plt.axis('equal')
                plt.title(section_name)
                plt.show()

        else:

            ops.section('Fiber', sec_tag)


            frt = -self.frc*(self.bf*self.tf) / \
                (self.bf*self.tf+self.tw*self.dw)
            # Define web fibers
            self.b = self.Hk/(self.E+self.Hk)
            ops.uniaxialMaterial('Steel01', matTag+1, self.Fy, self.E, self.b)
            ops.uniaxialMaterial('InitStressMaterial',
                                    matTag, matTag+1, frt)
            ops.patch('rect', matTag,Nfw,1,-self.tw/2,-self.dw/2, self.tw/2,self.dw/2)

            # Define flange fibers
            region_width = self.bf/self.num_regions
            Nff = ceil(region_width*(self.num_fiber/self.bf))
            for i in range(self.num_regions):
                fri = self.frc + ((i+0.5)/self.num_regions)*(frt-self.frc)

                matTagi = matTag+2*(i+1)
                ops.uniaxialMaterial(
                    'Steel01', matTagi+1, self.Fy, self.E, self.b)
                ops.uniaxialMaterial(
                    'InitStressMaterial', matTagi, matTagi+1, fri)

                ops.patch('rect', matTagi,Nff, 1,  -region_width/2, -self.d/2, region_width/2,-self.dw/2)
                ops.patch('rect', matTagi, Nff, 1,  -region_width/2, self.dw/2, region_width/2,self.d/2)
                
            if plot:

                fib_sec = [

                        ['section', 'Fiber', sec_tag, '-GJ', 1.0e4],

                        ["patch","rect", matTag,Nfw,1,-self.tw/2,-self.dw/2, self.tw/2,self.dw/2],
                        ["patch","rect",matTagi, Nff, 1,  -region_width/2, -self.d/2, region_width/2,-self.dw/2],
                        ["patch","rect", matTagi, Nff, 1,  -region_width/2, self.dw/2, region_width/2,self.d/2]

                        ]
                
                matcolor = [ 'gold']*10000
                opsv.plot_fiber_section(fib_sec, matcolor=matcolor)
                plt.axis('equal')
                plt.title(section_name)
                plt.show()
