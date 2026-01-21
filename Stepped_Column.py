import openseespy.opensees as ops
from libdenavit.section.wide_flange import *
from libdenavit.OpenSees import AnalysisResults
from Units import *  
import opsvis as opsv
from libdenavit.OpenSees.plotting import *
from libdenavit import find_limit_point_in_list, interpolate_list
from Plots import line_plot
import numpy as np

class Stepped_Column:

    def __init__(self,bottom_column_section_name,height_of_bottom_column,number_of_elements_bottom_column,load_on_bottom_column,
                 top_column_section_name,height_of_top_column,number_of_elements_top_column,load_on_top_column,offset_1,offset_3,**kwargs):

        self.height_of_bottom_column = height_of_bottom_column
        self.height_of_top_column = height_of_top_column
        self.number_of_elements_bottom_column = number_of_elements_bottom_column
        self.number_of_elements_top_column = number_of_elements_top_column
        self.load_on_bottom_column = load_on_bottom_column
        self.load_on_top_column = load_on_top_column
        self.bottom_column_section_name = bottom_column_section_name
        self.top_column_section_name = top_column_section_name 
        self.offset1 = offset_1  ## dist from top of bottom column to right node where point load is applied for bottom column
        self.offset3 = offset_3  ## dist from top of top column to right node where load is applied for top column
        self.length_of_single_bottom_column_element = height_of_bottom_column/number_of_elements_bottom_column 
        self.length_of_single_top_column_element = height_of_top_column/number_of_elements_top_column
        self.total_height = height_of_bottom_column + height_of_top_column

        defaults={'nip':3,
                  'mat_type':'Steel01',
                  'nfy':20,
                  'nfx':20,
                  'num_regions':10,
                  'Fy':36,
                  'E':29000,
                  'frc':0.0,
                  'Elastic_analysis':False,
                  'Second_order_effects':True,
                  'Residual_Stress':True,
                  'Geometric_Imperfection':False,
                  'geometric_imperfection_ratio':1/500,
                  'stiffness_reduction':1,
                  'strength_reduction':1,
                  'plot_model':True}    
        
        for key,value in defaults.items():
            setattr(self,key,kwargs.get(key,value))


    def build_stepped_column(self):

        ops.wipe()
        ops.model('basic','-ndm',2,'-ndf',3)

        ## Calculate left offset based on section depth of top and bottom column
        bottom_column_section=WF_Database(self.bottom_column_section_name,unit=1) 
        top_column_section=WF_Database(self.top_column_section_name,unit=1)
        self.offset2= (bottom_column_section.d/2)-(top_column_section.d/2)     
        print("Offset2 (Left offset for top column): ",self.offset2)   
        input()

        ## Bottom node of bottom column
        ops.node(1,0.0,0.0)

        ## Create nodes along bottom column
        for i in range(self.number_of_elements_bottom_column):
            x=0.0
            y=(i+1)*self.length_of_single_bottom_column_element
            x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
            ops.node(i+2,x + x_imp,y)
        self.bottom_column_top_node_tag=self.number_of_elements_bottom_column+1

        ## Create right offset node for bottom column
        right_offset1_node_tag=self.bottom_column_top_node_tag+1
        x=self.offset1
        y=self.height_of_bottom_column
        x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
        ops.node(right_offset1_node_tag,x + x_imp,y)


        ## Create nodes along top column
        self.top_column_bottom_node_tag=right_offset1_node_tag+1
        y = self.height_of_bottom_column  # bottom of top column
        x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
        ops.node(self.top_column_bottom_node_tag,-self.offset2 + x_imp,self.height_of_bottom_column)
        for j in range(self.number_of_elements_top_column):
            y = self.height_of_bottom_column+(j+1)*self.length_of_single_top_column_element
            x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
            ops.node(self.top_column_bottom_node_tag+j+1,-self.offset2 + x_imp, y)
        self.top_column_top_node_tag=self.top_column_bottom_node_tag+self.number_of_elements_top_column

        ## Create right offset node for top column
        right_offset3_node_tag=self.top_column_top_node_tag+1
        x=self.offset3-self.offset2
        y=self.height_of_bottom_column+self.height_of_top_column
        x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
        ops.node(right_offset3_node_tag,x + x_imp,y)
        

        ## Define geometric transformation for columns
        col_TransTag=1
        if self.Second_order_effects:    
            ops.geomTransf("PDelta", col_TransTag)

        else:
            ops.geomTransf("Linear", col_TransTag)

        if self.Elastic_analysis:
            mat_type = 'Elastic'
            frc = 0
        else:
            mat_type = self.mat_type
            frc = -0.3 * self.Fy if self.Residual_Stress else 0

        ## Fix base of bottom column
        ops.fix(1,1,1,0)

        ## Define bottom column section 
        bottom_column_section_tag=1
        
        bottom_column = I_shape(bottom_column_section.d, bottom_column_section.tw, bottom_column_section.bf, bottom_column_section.tf,   
            self.Fy,self.E,
            A=bottom_column_section.A, 
            Ix=bottom_column_section.Ix,Zx=bottom_column_section.Zx,Sx=bottom_column_section.Sx,rx=bottom_column_section.rx,
            Iy=bottom_column_section.Iy,Zy=bottom_column_section.Zy,Sy=bottom_column_section.Sy,ry=bottom_column_section.ry,
            J=bottom_column_section.J,Cw=bottom_column_section.Cw,rts=bottom_column_section.rts,ho=bottom_column_section.ho)
        bottom_column.build_ops_fiber_section(bottom_column_section_tag,
                                    start_material_id=1,
                                    mat_type=mat_type,
                                    nfy=self.nfy, nfx=self.nfx,
                                    frc=frc,num_regions=self.num_regions,
                                    stiffness_reduction=self.stiffness_reduction,strength_reduction=self.strength_reduction,
                                    axis='x')
        
        ## Define bottom column elements
        ops.beamIntegration("Lobatto", bottom_column_section_tag, bottom_column_section_tag, self.nip)
        for k in range(self.number_of_elements_bottom_column):
            element_tag=k+1
            node_i_tag=k+1
            node_j_tag=k+2
            ops.element('forceBeamColumn',element_tag,node_i_tag,node_j_tag,
                        col_TransTag,
                        bottom_column_section_tag,'-mass', 1)
        
        ## Define offset beam section
        offset_beam_section_tag=2
        ops.section('Elastic', offset_beam_section_tag, 29000*ksi, 1000*(inch**2), 1.0e6*(inch**4))
        ## Define beam integration for offset beam section
        ops.beamIntegration("Lobatto", offset_beam_section_tag, offset_beam_section_tag, self.nip)
        ## Define offset beam element at top of bottom column
        offset1_element_tag=self.number_of_elements_bottom_column+1
        ops.element('forceBeamColumn',offset1_element_tag,self.bottom_column_top_node_tag,right_offset1_node_tag,
                    col_TransTag,offset_beam_section_tag,'-mass', 1)
        
        ## Define 2nd offset beam element at bottom of top column
        offset2_element_tag=offset1_element_tag+1
        ops.element('forceBeamColumn',offset2_element_tag,self.bottom_column_top_node_tag,self.top_column_bottom_node_tag,
                    col_TransTag,offset_beam_section_tag,'-mass', 1)
        
        ## Define top column section
        top_column_section_tag=3

        top_column = I_shape(top_column_section.d, top_column_section.tw, top_column_section.bf, top_column_section.tf,   
            self.Fy,self.E,
            A=top_column_section.A, 
            Ix=top_column_section.Ix,Zx=top_column_section.Zx,Sx=top_column_section.Sx,rx=top_column_section.rx,
            Iy=top_column_section.Iy,Zy=top_column_section.Zy,Sy=top_column_section.Sy,ry=top_column_section.ry,
            J=top_column_section.J,Cw=top_column_section.Cw,rts=top_column_section.rts,ho=top_column_section.ho)
        top_column.build_ops_fiber_section(top_column_section_tag,
                                    start_material_id=200,
                                    mat_type=mat_type,
                                    nfy=self.nfy, nfx=self.nfx,
                                    frc=frc,num_regions=self.num_regions,
                                    stiffness_reduction=self.stiffness_reduction,strength_reduction=self.strength_reduction,
                                    axis='x')
        ## Define beam integration for top column section
        ops.beamIntegration("Lobatto", top_column_section_tag, top_column_section_tag, self.nip)       
        ## Define top column elements
        for m in range(self.number_of_elements_top_column):
            element_tag=self.number_of_elements_bottom_column+3+m
            node_i_tag=self.top_column_bottom_node_tag+m
            node_j_tag=self.top_column_bottom_node_tag+m+1
            ops.element('forceBeamColumn',element_tag,node_i_tag,node_j_tag,
                        col_TransTag,
                        top_column_section_tag,'-mass', 1)
            
        ## Define offset beam element at top of top column
        offset3_element_tag=self.number_of_elements_bottom_column+3+self.number_of_elements_top_column 
        ops.element('forceBeamColumn',offset3_element_tag,self.top_column_top_node_tag,right_offset3_node_tag,
                    col_TransTag,offset_beam_section_tag,'-mass', 1)
        
        ## Support at top node of top column 
        ops.fix(self.top_column_top_node_tag,1,0,0)

        ## Apply point load at right offset node of bottom column and offset node of top column
        ops.timeSeries('Linear',1)
        ops.pattern('Plain',1,1)
        ops.load(right_offset1_node_tag,0.0,-self.load_on_bottom_column,0.0)  ## Apply vertical downward load
        ops.load(right_offset3_node_tag,0.0,-self.load_on_top_column,0.0)  ## Apply vertical downward load

        if self.plot_model: 
            opsv.plot_model()
            opsv.plot_load()
    
    
    def show_model(self):
        # show quick model diagnostics then plot using libdenavit's plotting
        try:
            node_coords = get_node_coords()
            element_nodes = get_element_nodes()
            print(f"Plotting model: {len(node_coords)} nodes, {len(element_nodes)} elements")
        except Exception:
            print("Unable to read node/element info for diagnostics")
        plot_undeformed_2d(axis_equal=True)

    def run_limit_point_anlaysis(self, target_disp=1, steps=10,**kwargs):
        incr_LCA= kwargs.get('incr_LCA', 0.1)          ######### LCA refers to Load Controlled Analysis
        num_steps_LCA= kwargs.get('num_steps_LCA', 10)            ######### LCA refers to Load Controlled Analysis
        steel_strain_limit = kwargs.get('steel_strain_limit', 0.05)
        eigenvalue_limit = kwargs.get('eigenvalue_limit', 0)
        try_smaller_steps = kwargs.get('try_smaller_steps', True)
        print_ops_status = kwargs.get('print_ops_status', True)
        # Initialize analysis results
        results = AnalysisResults()
        attributes = ['load_ratio','control_node_displacement',
                      'lowest_eigenvalue']
        
        for attr in attributes:
            setattr(results, attr, [])

        # Define function to find limit point
        def find_limit_point():
            if 'Analysis Failed' in results.exit_message:
                print('Analysis Failed before reaching negative eigenvalue')
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif 'Eigenvalue Limit' in results.exit_message:
                ind, x = find_limit_point_in_list(results.lowest_eigenvalue, eigenvalue_limit)
            elif 'Extreme Steel Fiber Strain Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.absolute_maximum_strain, steel_strain_limit)
            elif  'Analysis Failed In Load Controlled Loading before entering Displacement controlled Loading' in results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))  
            elif 'Moving to Displacement Controlled Analysis' in results.exit_message:
                print('Moving to Displacement Controlled Analysis, no limit point reached in Load Controlled Analysis')
                return          
            else:
                raise Exception('Unknown limit point')
            results.maximum_load_ratio_at_limit_point = interpolate_list(results.load_ratio, ind, x)
            print(' Max Load Ratio',results.maximum_load_ratio_at_limit_point)

        control_node = self.bottom_column_top_node_tag  # Bottom node of top column is the control node for displacement control
        # print(f'Control node for displacement control: {control_node}')
        control_dof = 1   # Horizontal displacement is the control DOF
        # region Define recorder
        def record():
            time = ops.getTime()
            results.load_ratio.append(time)
            results.lowest_eigenvalue.append(ops.eigen("-genBandArpack", 1)[0])
            results.control_node_displacement.append(ops.nodeDisp(control_node, control_dof))
        # endregion

        ops.initialize()
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 10)
        ops.algorithm('Newton')
    
        ops.integrator('LoadControl', incr_LCA)  
        ops.analysis('Static')
        record()
        for i in range(num_steps_LCA):
            if print_ops_status:
                print(f'Running Load Controlled Analysis Step {i}')
            ok = ops.analyze(1)
            if ok != 0:
                print(f'Load controlled analysis failed in step {i}')
                results.exit_message = 'Analysis Failed In Load Controlled Loading before entering Displacement controlled Loading'
                find_limit_point()
                return results
            else:
                print('Load controlled analysis PASSED')
                results.exit_message='Moving to Displacement Controlled Analysis'
            record()
        
        print('Starting Displacement Controlled Analysis')
        ops.loadConst('-time', 0.0)
        ops.timeSeries('Linear', 2)
        ops.pattern('Plain',2,2)
        ## Apply lateral load at the control node
        ops.load(control_node,1000,0,0)  ## Apply horizontal load 
        dU = target_disp / steps
        ops.integrator('DisplacementControl', control_node, control_dof, dU)

        record()
        i=1
        while True:
            print(f'Running Displacement Controlled Analysis {i}')
            i=i+1
            ok = ops.analyze(1)
            if try_smaller_steps:
                if ok != 0:
                    if print_ops_status:
                        print(f'Trying the step size of: {dU / 10}')
                    ops.integrator('DisplacementControl',control_node, control_dof, dU / 10)
                    ok = ops.analyze(1)

                if ok != 0:
                    if print_ops_status:
                        print(f'Trying the step size of: {dU / 100}')
                    ops.integrator('DisplacementControl',control_node, control_dof, dU / 100)
                    ok = ops.analyze(1)

                if ok != 0:
                    if print_ops_status:
                        print(f'Trying the step size of: {dU / 1000}')
                    ops.integrator('DisplacementControl', control_node, control_dof, dU / 1000)
                    ok = ops.analyze(1)
                    if ok == 0:
                        dU = dU / 10
                        if print_ops_status:
                            print(f'Changed the step size to: {dU}')

                if ok != 0:
                    if print_ops_status:
                        print(f'Trying the step size of: {dU / 10000}')
                    ops.integrator('DisplacementControl', control_node, control_dof, dU / 10000)
                    ok = ops.analyze(1)
                    if ok == 0:
                        dU = dU / 10
                        if print_ops_status:
                            print(f'Changed the step size to: {dU / 10}')

            if ok != 0:
                if print_ops_status:
                    print('Trying ModifiedNewton')
                ops.algorithm('ModifiedNewton')
                ok = ops.analyze(1)
                if ok == 0:
                    if print_ops_status:
                        print('ModifiedNewton worked')

            if ok != 0:
                if print_ops_status:
                    print('Trying KrylovNewton')
                ops.algorithm('KrylovNewton')
                ok = ops.analyze(1)
                if ok == 0:
                    if print_ops_status:
                        print('KrylovNewton worked')

            if ok != 0:
                if print_ops_status:
                    print('Trying KrylovNewton and Greater Tolerance')
                ops.algorithm('KrylovNewton')
                ops.test('NormUnbalance', 1e-4, 10)
                ok = ops.analyze(1)
                if ok == 0:
                    if print_ops_status:
                        print('KrylovNewton worked')

            if ok == 0:
                # Reset analysis options
                ops.algorithm('Newton')
                ops.test('NormUnbalance', 1e-3, 10)
                ops.integrator('DisplacementControl', control_node, control_dof, dU)
            else:
                print('Analysis Failed')
                results.exit_message = 'Analysis Failed'
                break


            record()

            # Check for lowest eigenvalue less than zero
            if eigenvalue_limit is not None:
                if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                    results.exit_message = 'Eigenvalue Limit Reached'
                    break

        find_limit_point()
        return results

Stepped_Column = Stepped_Column(bottom_column_section_name='W14X132',height_of_bottom_column=4.0*ft,number_of_elements_bottom_column=5,load_on_bottom_column=100*kip,
                                    top_column_section_name='W14X90',height_of_top_column=3.0*ft,number_of_elements_top_column=5,load_on_top_column=100*kip,
                                    offset_1=0.3,offset_3=0.25,
                                    Fy=36*ksi,E=29000*ksi,
                                    Elastic_analysis=False,
                                    Second_order_effects=True,
                                    Residual_Stress=True,
                                    Geometric_Imperfection=True,
                                    geometric_imperfection_ratio=1/10)

Stepped_Column.build_stepped_column()
Stepped_Column.show_model()
results = Stepped_Column.run_limit_point_anlaysis(target_disp=10, steps=100)

line_plot(results.control_node_displacement, results.load_ratio,
        xlabel='Displacement at Control Node', ylabel='Load Ratio λ',
        title='Load Ratio vs Displacement',show=True)

line_plot(results.lowest_eigenvalue, results.load_ratio,
        xlabel='Lowest Eigenvalue', ylabel='Load Ratio λ',
        title='Load Ratio vs Lowest Eigenvalue',show=True)