import openseespy.opensees as ops
from libdenavit.section.wide_flange import I_shape
from libdenavit.OpenSees import AnalysisResults
import Units as U
import opsvis as opsv
from libdenavit.OpenSees.plotting import plot_undeformed_2d,get_element_nodes,get_node_coords
from libdenavit import find_limit_point_in_list, interpolate_list
from Plots import line_plot
import numpy as np
from Structures_2D import Structures_2D

class Stepped_Column(Structures_2D):

    def __init__(self,bottom_column_section_name,height_of_bottom_column,load_on_bottom_column,
                 top_column_section_name,height_of_top_column,number_of_elements,load_on_top_column,lateral_load,offset_1,offset_3,**kwargs):

        self.height_of_bottom_column = height_of_bottom_column
        self.height_of_top_column = height_of_top_column
        self.no_of_elements_column = number_of_elements
        self.load_on_bottom_column = load_on_bottom_column
        self.load_on_top_column = load_on_top_column
        self.lateral_load=lateral_load
        self.bottom_column_section_name = bottom_column_section_name
        self.top_column_section_name = top_column_section_name 
        self.offset1 = offset_1  ## dist from top of bottom column to right node where point load is applied for bottom column
        self.offset3 = offset_3  ## dist from top of top column to right node where load is applied for top column
        self.length_of_single_bottom_column_element = height_of_bottom_column/number_of_elements
        self.length_of_single_top_column_element = height_of_top_column/number_of_elements
        self.total_height = height_of_bottom_column + height_of_top_column
        self.all_element_connectivity_section_and_bending_axes_detail = []
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
                  'plot_model':True,
                  'bottom_column_bending_axes':'x',
                  'top_column_bending_axes':'x'
                  }    
        
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
        for i in range(self.no_of_elements_column):
            x=0.0
            y=(i+1)*self.length_of_single_bottom_column_element
            x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
            ops.node(i+2,x + x_imp,y)
        self.bottom_column_top_node_tag=self.no_of_elements_column+1

        ## Create right offset node for bottom column
        self.right_offset1_node_tag=self.bottom_column_top_node_tag+1
        x=self.offset1
        y=self.height_of_bottom_column
        x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
        ops.node(self.right_offset1_node_tag,x + x_imp,y)


        ## Create nodes along top column
        self.top_column_bottom_node_tag=self.right_offset1_node_tag+1
        y = self.height_of_bottom_column  # bottom of top column
        x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
        ops.node(self.top_column_bottom_node_tag,-self.offset2 + x_imp,self.height_of_bottom_column)
        for j in range(self.no_of_elements_column):
            y = self.height_of_bottom_column+(j+1)*self.length_of_single_top_column_element
            x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
            ops.node(self.top_column_bottom_node_tag+j+1,-self.offset2 + x_imp, y)
        self.top_column_top_node_tag=self.top_column_bottom_node_tag+self.no_of_elements_column

        ## Create right offset node for top column
        self.right_offset3_node_tag=self.top_column_top_node_tag+1
        x=self.offset3-self.offset2
        y=self.height_of_bottom_column+self.height_of_top_column
        x_imp = self.geometric_imperfection_ratio * np.sin(np.pi * y / self.total_height)
        ops.node(self.right_offset3_node_tag,x + x_imp,y)
        

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
        setattr(self,self.bottom_column_section_name,bottom_column)
        ## Define bottom column elements
        ops.beamIntegration("Lobatto", bottom_column_section_tag, bottom_column_section_tag, self.nip)
        for k in range(self.no_of_elements_column):
            element_tag=k+1
            node_i_tag=k+1
            node_j_tag=k+2
            ops.element('forceBeamColumn',element_tag,node_i_tag,node_j_tag,
                        col_TransTag,
                        bottom_column_section_tag,'-mass', 1)
            self.all_element_connectivity_section_and_bending_axes_detail.append([element_tag,node_i_tag,node_j_tag,self.bottom_column_section_name,self.bottom_column_bending_axes,'col'])
        
        ## Define offset beam section
        offset_beam_section_tag=2
        ops.section('Elastic', offset_beam_section_tag, 29000*ksi, 1000*(inch**2), 1.0e6*(inch**4))
        ## Define beam integration for offset beam section
        ops.beamIntegration("Lobatto", offset_beam_section_tag, offset_beam_section_tag, self.nip)
        ## Define offset beam element at top of bottom column
        offset1_element_tag=self.no_of_elements_column+1
        ops.element('forceBeamColumn',offset1_element_tag,self.bottom_column_top_node_tag,self.right_offset1_node_tag,
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
        setattr(self,self.top_column_section_name,top_column)
        ## Define beam integration for top column section
        ops.beamIntegration("Lobatto", top_column_section_tag, top_column_section_tag, self.nip)       
        ## Define top column elements
        for m in range(self.no_of_elements_column):
            element_tag=self.no_of_elements_column+3+m
            node_i_tag=self.top_column_bottom_node_tag+m
            node_j_tag=self.top_column_bottom_node_tag+m+1
            ops.element('forceBeamColumn',element_tag,node_i_tag,node_j_tag,
                        col_TransTag,
                        top_column_section_tag,'-mass', 1)
            self.all_element_connectivity_section_and_bending_axes_detail.append([element_tag,node_i_tag,node_j_tag,self.top_column_section_name,self.top_column_bending_axes,'col'])
        ## Define offset beam element at top of top column
        offset3_element_tag=self.no_of_elements_column+3+self.no_of_elements_column 
        ops.element('forceBeamColumn',offset3_element_tag,self.top_column_top_node_tag,self.right_offset3_node_tag,
                    col_TransTag,offset_beam_section_tag,'-mass', 1)
        
        ## Support at top node of top column 
        ops.fix(self.top_column_top_node_tag,1,0,0)


        if self.plot_model: 
            opsv.plot_model()
            opsv.plot_load()

    def add_vertical_load(self) :
        ops.load(self.right_offset1_node_tag,0.0,-self.load_on_bottom_column,0.0)  ## Apply vertical downward load
        ops.load(self.right_offset3_node_tag,0.0,-self.load_on_top_column,0.0)  ## Apply vertical downward load

    def add_lateral_load(self):
        ops.load(self.bottom_column_top_node_tag,self.lateral_load,0,0)  ## Apply horizontal load

    def show_model(self):
        # show quick model diagnostics then plot using libdenavit's plotting
        try:
            node_coords = get_node_coords()
            element_nodes = get_element_nodes()
            print(f"Plotting model: {len(node_coords)} nodes, {len(element_nodes)} elements")
        except Exception:
            print("Unable to read node/element info for diagnostics")
        plot_undeformed_2d(axis_equal=True)

    def run_load_controlled_analysis(self, **kwargs):
        incr_LCA= kwargs.get('incr_LCA', 0.1)          ######### LCA refers to Load Controlled Analysis
        num_steps_LCA= kwargs.get('num_steps_LCA', 10)            ######### LCA refers to Load Controlled Analysis
        steel_strain_limit = kwargs.get('steel_strain_limit', 0.05)
        eigenvalue_limit = kwargs.get('eigenvalue_limit', 0)
        P_M_M_interaction_limit=kwargs.get('P_M_M_interaction_limit',1)
        try_smaller_steps = kwargs.get('try_smaller_steps', True)
        print_ops_status = kwargs.get('print_ops_status', True)
        # Initialize analysis results
        results = AnalysisResults()
        attributes = ['load_ratio','control_node_displacement',
                      'lowest_eigenvalue','vertical_reaction','lateral_reaction','absolute_maximum_strain','max_P_M_M_interaction']
        
        for attr in attributes:
            setattr(results, attr, [])

        # Define function to find limit point
        def find_limit_point():
            if 'Analysis Failed' in results.exit_message:
                print('Analysis Failed before reaching negative eigenvalue')
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif 'Eigenvalue Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.lowest_eigenvalue, eigenvalue_limit)
            elif 'Extreme Steel Fiber Strain Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.absolute_maximum_strain, steel_strain_limit)
            elif  'Analysis Failed In Load Controlled Loading before entering Displacement controlled Loading' in results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))  
            elif 'Full Load Applied' in results.exit_message:
                print('Full Load Applied')
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif 'P_M_M interaction Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.max_P_M_M_interaction, P_M_M_interaction_limit)          
            else:
                raise Exception('Unknown limit point')
            results.maximum_load_ratio_at_limit_point = interpolate_list(results.load_ratio, ind, x)
            print(' Max Load Ratio',results.maximum_load_ratio_at_limit_point)


        control_node = self.bottom_column_top_node_tag  # Bottom node of top column is the control node for displacement control
        # print(f'Control node for displacement control: {control_node}')
        control_dof = 1   # Horizontal displacement is the control DOF

        def record():
            time = ops.getTime()
            results.load_ratio.append(time)
            results.lowest_eigenvalue.append(ops.eigen("-genBandArpack", 1)[0])
            results.control_node_displacement.append(ops.nodeDisp(control_node, control_dof))
            ops.reactions()
            total_vertical_rxn=ops.nodeReaction(1)[1] 
            lateral_reaction=ops.nodeReaction(1)[0] + ops.nodeReaction(self.top_column_top_node_tag)[0]
            results.vertical_reaction.append(total_vertical_rxn)
            results.lateral_reaction.append(lateral_reaction)            
            results.absolute_maximum_strain.append(self.return_max_of_fiber_strain_in_all_elements())
            # results.control_node_displacement.append(ops.nodeDisp(control_node, control_dof))
            max_PMM, max_ele_tag,P_M_M_interaction_all_elements,Element_Forces=self.return_P_M_M_interaction_values()
            results.max_P_M_M_interaction.append(max_PMM)

        ops.initialize()
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 10)
        ops.algorithm('Newton')
    
        ops.timeSeries('Linear',1)
        ops.pattern('Plain',1,1)
        self.add_vertical_load()
        self.add_lateral_load()

        ops.integrator('LoadControl', 1/num_steps_LCA)  
        ops.analysis('Static')
        record()
        for i in range(num_steps_LCA):
            if print_ops_status:
                print(f'Running Load Controlled Analysis Step {i}')
            ok = ops.analyze(1)
            if ok != 0:
                print(f'Load controlled analysis failed in step {i}')
                results.exit_message = 'Analysis Failed In Load Controlled Loading '
                find_limit_point()
                
            else:
                print('Load controlled analysis PASSED')
                results.exit_message='Full Load Applied.'
                
            record()
            # Check for lowest eigenvalue less than zero
            if eigenvalue_limit is not None:
                if results.lowest_eigenvalue[-1] < eigenvalue_limit:
                    results.exit_message = 'Eigenvalue Limit Reached'
                    find_limit_point()
                    return results
                    # break

            # Check for strain in extreme steel fiber
            if steel_strain_limit is not None:
                # if Structures_2D.print_ops_status:
                #     print(f'Checking Steel Tensile Strain')
                if results.absolute_maximum_strain[-1] > steel_strain_limit:
                    results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                    find_limit_point()
                    return results
                    # break
            # Check for maximum PMM interaction value    
            if self.Elastic_analysis:
                if P_M_M_interaction_limit is not None:
                    # if Structures_2D.print_ops_status:
                    #     print(f'Checking PMM Interaction')
                    if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                        results.exit_message = 'P_M_M interaction Limit Reached'
                        find_limit_point()
                        return results
        find_limit_point()
        return results



    def run_displacement_controlled_analysis(self, target_disp=1, steps=10,**kwargs):
        incr_LCA= kwargs.get('incr_LCA', 0.1)          ######### LCA refers to Load Controlled Analysis
        num_steps_LCA= kwargs.get('num_steps_LCA', 10)            ######### LCA refers to Load Controlled Analysis
        steel_strain_limit = kwargs.get('steel_strain_limit', 0.05)
        eigenvalue_limit = kwargs.get('eigenvalue_limit', 0)
        P_M_M_interaction_limit=kwargs.get('P_M_M_interaction_limit',1)
        try_smaller_steps = kwargs.get('try_smaller_steps', True)
        print_ops_status = kwargs.get('print_ops_status', True)
        # Initialize analysis results
        results = AnalysisResults()
        attributes =['load_ratio','control_node_displacement',
                      'lowest_eigenvalue','vertical_reaction','lateral_reaction','absolute_maximum_strain','max_P_M_M_interaction']
        
        for attr in attributes:
            setattr(results, attr, [])

        # Define function to find limit point
        def find_limit_point():
            print(results.exit_message)
            if 'Analysis Failed' in results.exit_message:
                print('Analysis Failed before reaching negative eigenvalue')
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif 'Eigenvalue Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.lowest_eigenvalue, eigenvalue_limit)
            elif 'Extreme Steel Fiber Strain Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.absolute_maximum_strain, steel_strain_limit)
            elif  'Analysis Failed In Load Controlled Loading before entering Displacement controlled Loading' in results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio)) 
            elif 'P_M_M interaction Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.max_P_M_M_interaction, P_M_M_interaction_limit)           
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
            ops.reactions()
            total_vertical_rxn=ops.nodeReaction(1)[1] 
            lateral_reaction=ops.nodeReaction(1)[0] + ops.nodeReaction(self.top_column_top_node_tag)[0]
            results.vertical_reaction.append(total_vertical_rxn)
            results.lateral_reaction.append(lateral_reaction)            
            results.absolute_maximum_strain.append(self.return_max_of_fiber_strain_in_all_elements())
            # results.control_node_displacement.append(ops.nodeDisp(control_node, control_dof))
            max_PMM, max_ele_tag,P_M_M_interaction_all_elements,Element_Forces=self.return_P_M_M_interaction_values()
            results.max_P_M_M_interaction.append(max_PMM)
        # endregion

        ops.initialize()
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 10)
        ops.algorithm('Newton')
    
        # ops.integrator('LoadControl', 1/num_steps_LCA)  
        ops.analysis('Static')
        record()
        # for i in range(num_steps_LCA):
        #     if print_ops_status:
        #         print(f'Running Load Controlled Analysis Step {i}')
        #     ok = ops.analyze(1)
        #     if ok != 0:
        #         print(f'Load controlled analysis failed in step {i}')
        #         results.exit_message = 'Analysis Failed In Load Controlled Loading before entering Displacement controlled Loading'
        #         find_limit_point()
        #         return results
        #     else:
        #         print('Load controlled analysis PASSED')
        #         results.exit_message='Moving to Displacement Controlled Analysis'
        #     record()
        
        print('Starting Displacement Controlled Analysis')
        # ops.loadConst('-time', 0.0)
        ops.timeSeries('Linear', 2)
        ops.pattern('Plain',2,2)
        self.add_vertical_load()
        self.add_lateral_load()
        ## Apply lateral load at the control node
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

            # Check for strain in extreme steel fiber
            if steel_strain_limit is not None:
                # if Structures_2D.print_ops_status:
                #     print(f'Checking Steel Tensile Strain')
                print(results.absolute_maximum_strain[-1])
                if results.absolute_maximum_strain[-1] > steel_strain_limit:
                    results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                    break
            # Check for maximum PMM interaction value    
            if self.Elastic_analysis:
                if P_M_M_interaction_limit is not None:
                    # if Structures_2D.print_ops_status:
                    #     print(f'Checking PMM Interaction')
                    if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                        results.exit_message = 'P_M_M interaction Limit Reached'
                        break
        find_limit_point()
        return results

if __name__ == "__main__":
    Stepped_Column = Stepped_Column(bottom_column_section_name='W14X132',height_of_bottom_column=8.0*U.ft,load_on_bottom_column=100*10*10*100*5*U.KN,
                                        top_column_section_name='W14X90',height_of_top_column=8.0*U.ft,number_of_elements=2,load_on_top_column=80*10*10*100*5*U.KN,lateral_load=20*10*10*100*5*U.KN,
                                        offset_1=0.3,offset_3=0.25,
                                        Fy=36*U.ksi,E=29000*U.ksi,
                                        Elastic_analysis=False,
                                        Second_order_effects=True,
                                        Residual_Stress=True,
                                        Geometric_Imperfection=True,
                                        geometric_imperfection_ratio=1/500)

    Stepped_Column.build_stepped_column()
    print(Stepped_Column.all_element_connectivity_section_and_bending_axes_detail)
    Stepped_Column.show_model()
    results= Stepped_Column.run_load_controlled_analysis()
    # results = Stepped_Column.run_displacement_controlled_analysis(target_disp=10, steps=100)

    line_plot( results.load_ratio,results.control_node_displacement,
            xlabel='Load Ratio λ', ylabel='Displacement at Control Node',
            title='Load Ratio vs Displacement',show=True)

    line_plot( results.load_ratio,results.lowest_eigenvalue,
            xlabel='Load Ratio λ', ylabel='Lowest Eigenvalue',
            title='Load Ratio vs Lowest Eigenvalue',show=True)

    line_plot( results.load_ratio,results.vertical_reaction, 
            xlabel='Load Ratio λ', ylabel='Vertical Reaction',
            title='Load Ratio vs Vertical Reaction',show=True)

    line_plot( results.load_ratio,results.lateral_reaction, 
            xlabel='Load Ratio λ', ylabel='Lateral Reaction',
            title='Load Ratio vs Lateral Reaction',show=True) 

    line_plot( results.load_ratio,results.absolute_maximum_strain, 
            xlabel='Load Ratio λ', ylabel='Absolute Maximum Strain',
            title='Load Ratio vs Absolute Maximum Strain',show=True) 

    line_plot( results.load_ratio,results.max_P_M_M_interaction,
            xlabel='Load Ratio λ', ylabel='Max P-M-M Interaction',
            title='Load Ratio vs Max P-M-M Interaction',show=True) 