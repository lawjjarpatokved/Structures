import math
import openseespy.opensees as ops
import opsvis as opsv
import matplotlib.pyplot as plt
import numpy as np

from libdenavit import find_limit_point_in_list, interpolate_list
from libdenavit.section import I_shape
from libdenavit.section.wide_flange import WideFlangeMember_AISC2022
from libdenavit.section.database import wide_flange_database
from libdenavit.OpenSees import AnalysisResults,plot_deformed_2d,plot_undeformed_2d

plt.set_loglevel("warning")


def return_strength_ratio(ele_dict):
    ### The ele_dict contains the details like the section object, bending axes, Lcx, Lcy, Lb etc for each element. We can use this to compute the strength ratio for each element and return the maximum strength ratio among all elements.
    ele_output = dict()
    max_PMM = 0
    max_ele_tag = ''
    
    for ele_tag,value in ele_dict.items():
        section_obj = value['section']
        ele_bending_axis = value['bending_axes']
        Lcx = value['Lcx']
        Lcy = value['Lcy']
        Lb = value['Lb']
        L = value['L']

        member_obj = WideFlangeMember_AISC2022(
            section_obj,
            Fy=section_obj.Fy,
            E=section_obj.E,
            L=L
        )

        # Required Strengths
        forces = ops.eleResponse(ele_tag, 'localForce')
        Pr = abs(forces[0])
        Mr_i = forces[2]
        Mr_j = forces[5]
        max_Mr = max(abs(Mr_i), abs(Mr_j))
        if ele_bending_axis == 'x':
            Mrx = max_Mr
            Mry = 0
        elif ele_bending_axis == 'y':
            Mrx = 0
            Mry = max_Mr
        else:
            raise ValueError("The axis is not supported.")

        # Available Strengths
        Pcc = member_obj.Pnc(Lcx=Lcx, Lcy=Lcy)
        Cb = 1 # @todo - For the future: update this somehow            
        Mcx = member_obj.Mnx(Lb=Lb, Cb=Cb)
        Mcy = member_obj.Mny()

        # Interaction Eqn H1-1a or H1-1b
        # @todo - What is the element is in tension?
        if Pr / Pcc >= 0.2:
            P_M_M_interaction = Pr / Pcc + (8 / 9) * ((Mrx / Mcx) + (Mry / Mcy))
        else:
            P_M_M_interaction = Pr / (2 * Pcc) + ((Mrx / Mcx) + (Mry / Mcy))

        # Update maximum
        if P_M_M_interaction > max_PMM:
            max_PMM = P_M_M_interaction
            max_ele_tag = ele_tag

        # Store element output
        ele_output[ele_tag] = {'P_M_M_interaction':P_M_M_interaction,
                               'Element_Forces':forces,
                              }

    return max_PMM, max_ele_tag, ele_output

def return_max_of_fiber_strain_in_all_elements(ele_dict):
    ## returns the maximum strain from among all elements in the Frame
    maximum_compression_strain=[]
    maximum_tensile_strain=[]
    for key,value in ele_dict.items():
        ele_tag = key
        section_obj = value['section']
        ele_bending_axis = value['bending_axes']

        compression_strain = []
        tensile_strain = []

        nip=len(ops.sectionLocation(ele_tag))
        for i in range(nip):
            axial_strain, curvatureX, curvatureY = 0, 0, 0
            if ele_bending_axis=='x':
                axial_strain, curvatureX = ops.eleResponse(ele_tag,  # element tag
                                                            'section', i+1,  # select integration point
                                                            'deformation')  # response type               
            elif ele_bending_axis=='y':
                axial_strain, curvatureY = ops.eleResponse(ele_tag,  # element tag
                                                            'section', i+1,  # select integration point
                                                            'deformation')  # response type
            else:
                raise ValueError("The axis is not supported.")
            
            compression_strain.append(section_obj.maximum_compression_strain(axial_strain,curvatureX,curvatureY))
            tensile_strain.append(section_obj.maximum_tensile_strain(axial_strain,curvatureX,curvatureY))
        maximum_compression_strain.append(max(compression_strain,key=lambda x:abs(x)))
        maximum_tensile_strain.append(max(tensile_strain,key=lambda x:abs(x)))
    max_abs_compression_strain= max(abs(c) for c in maximum_compression_strain)
    max_abs_tensile_strain= max(abs(t) for t in maximum_tensile_strain)
    return max(max_abs_compression_strain,max_abs_tensile_strain)

class Stepped_Column():

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
        self.element_dict={}

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
                  'plot_model':False,
                  'bottom_column_bending_axes':'x',
                  'top_column_bending_axes':'x'
                  }    
        
        for key,value in defaults.items():
            setattr(self,key,kwargs.get(key,value))


    @property
    def total_height(self):
        return self.height_of_bottom_column + self.height_of_top_column

    @property
    def offset2(self):
        # Calculate left offset based on section depth of top and bottom column
        bottom_column_section_data = wide_flange_database[self.bottom_column_section_name]
        top_column_section_data = wide_flange_database[self.top_column_section_name]
        return (bottom_column_section_data['d']/2)-(top_column_section_data['d']/2)  
        
    def build_stepped_column(self):

        # Build model
        ops.wipe()
        ops.model('basic','-ndm',2,'-ndf',3)

        ## Bottom node of bottom column
        ops.node(1,0.0,0.0)

        ## Create nodes along bottom column
        for i in range(self.no_of_elements_column):
            x=0.0
            y=(i+1)*self.height_of_bottom_column/self.no_of_elements_column
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
            y = self.height_of_bottom_column+(j+1)*self.height_of_top_column/self.no_of_elements_column
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
            ops.geomTransf("Corotational", col_TransTag)

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
        bottom_column = I_shape.from_database(self.bottom_column_section_name,self.Fy,self.E)
        bottom_column.build_ops_fiber_section(bottom_column_section_tag,
                                    start_material_id=1,
                                    mat_type=mat_type,
                                    nfy=self.nfy, nfx=self.nfx,
                                    frc=frc,num_regions=self.num_regions,
                                    stiffness_reduction=self.stiffness_reduction,strength_reduction=self.strength_reduction,
                                    axis='x')

        
        ## Define bottom column elements
        ops.beamIntegration("Lobatto", bottom_column_section_tag, bottom_column_section_tag, self.nip)
        for k in range(self.no_of_elements_column):
            element_tag=k+1
            node_i_tag=k+1
            node_j_tag=k+2
            ops.element('forceBeamColumn',element_tag,node_i_tag,node_j_tag,
                        col_TransTag,
                        bottom_column_section_tag,'-mass', 1)
            coords_i = ops.nodeCoord(node_i_tag)
            coords_j = ops.nodeCoord(node_j_tag)
            L = math.sqrt((coords_j[0] - coords_i[0])**2 + (coords_j[1] - coords_i[1])**2)
            Leff = L * (self.no_of_elements_column)
            if self.bottom_column_bending_axes == 'x':
                Lcx = Leff
                Lcy = 0
            elif self.bottom_column_bending_axes == 'y':
                Lcx = 0
                Lcy = Leff
            Lb=0
            self.element_dict[element_tag] = {'section':bottom_column,'bending_axes':self.bottom_column_bending_axes,'Lcx':Lcx,'Lcy':Lcy,'Lb':Lb,'L':L}

        ## Define offset beam section
        offset_beam_section_tag=2
        ops.section('Elastic', offset_beam_section_tag, 29000, 1000, 1.0e6)
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
        top_column = I_shape.from_database(self.top_column_section_name,self.Fy,self.E)
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
        for m in range(self.no_of_elements_column):
            element_tag=self.no_of_elements_column+3+m
            node_i_tag=self.top_column_bottom_node_tag+m
            node_j_tag=self.top_column_bottom_node_tag+m+1
            ops.element('forceBeamColumn',element_tag,node_i_tag,node_j_tag,
                        col_TransTag,
                        top_column_section_tag,'-mass', 1)
            coords_i = ops.nodeCoord(node_i_tag)
            coords_j = ops.nodeCoord(node_j_tag)
            L = math.sqrt((coords_j[0] - coords_i[0])**2 + (coords_j[1] - coords_i[1])**2)
            Leff = L * (self.no_of_elements_column)
            if self.top_column_bending_axes == 'x':
                Lcx = Leff
                Lcy = 0
            elif self.top_column_bending_axes == 'y':
                Lcx = 0
                Lcy = Leff
            Lb=0
            self.element_dict[element_tag] = {'section':top_column,'bending_axes':self.top_column_bending_axes,'Lcx':Lcx,'Lcy':Lcy,'Lb':Lb,'L':L}  
               
        ## Define offset beam element at top of top column
        offset3_element_tag=self.no_of_elements_column+3+self.no_of_elements_column 
        ops.element('forceBeamColumn',offset3_element_tag,self.top_column_top_node_tag,self.right_offset3_node_tag,
                    col_TransTag,offset_beam_section_tag,'-mass', 1)
        
        ## Support at top node of top column 
        ops.fix(self.top_column_top_node_tag,1,0,0)


        # Define loads
        ops.timeSeries('Linear',1)
        ops.pattern('Plain',1,1)
        ops.load(self.right_offset1_node_tag,0.0,-self.load_on_bottom_column,0.0)
        ops.load(self.right_offset3_node_tag,0.0,-self.load_on_top_column,0.0)
        ops.load(self.bottom_column_top_node_tag,self.lateral_load,0,0)

        # Display model if requested
        if self.plot_model: 
            opsv.plot_model()
            opsv.plot_load()

    def show_model(self):
        plot_undeformed_2d(axis_equal=True)

    def run_load_controlled_analysis(self, target_load_factor=1.0, steps=1000, **kwargs):
        steel_strain_limit = kwargs.get('steel_strain_limit', None)
        eigenvalue_limit = kwargs.get('eigenvalue_limit', None)
        P_M_M_interaction_limit=kwargs.get('P_M_M_interaction_limit',1)
        print_ops_status = kwargs.get('print_ops_status', True)

        # Define control node
        control_node = self.bottom_column_top_node_tag  # Bottom node of top column is the control node for displacement control
        print(f'Control node for displacement control: {control_node}')
        control_dof = 1   # Horizontal displacement is the control DOF

        # Initialize analysis results
        results = AnalysisResults(initialize_empty_lists = ['load_ratio',
                                                            'control_node_displacement',
                                                            'lowest_eigenvalue',
                                                            'vertical_reaction',
                                                            'lateral_reaction',
                                                            'absolute_maximum_strain',
                                                            'max_P_M_M_interaction',
                                                            ])

        # Define function to find limit point
        def find_limit_point():
            if 'Analysis Failed' in results.exit_message:
                print('Analysis Failed before reaching negative eigenvalue')
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio))
            elif 'Eigenvalue Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.lowest_eigenvalue, eigenvalue_limit)
            elif 'Extreme Steel Fiber Strain Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.absolute_maximum_strain, steel_strain_limit)
            elif 'P_M_M interaction Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.max_P_M_M_interaction, P_M_M_interaction_limit)          
            elif 'Full Load Applied' in results.exit_message:
                results.maximum_load_ratio_at_limit_point = []
                return
            else:
                raise Exception('Unknown limit point')
            
            results.maximum_load_ratio_at_limit_point = interpolate_list(results.load_ratio, ind, x)
            print(' Max Load Ratio',results.maximum_load_ratio_at_limit_point)

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
            results.absolute_maximum_strain.append(return_max_of_fiber_strain_in_all_elements(self.element_dict))
            max_PMM, max_ele_tag,_ = return_strength_ratio(self.element_dict)
            results.max_P_M_M_interaction.append(max_PMM)

        # Define analysis options
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 10, 1)
        ops.algorithm('Newton')   
        ops.analysis('Static')

        # Run one step with no load
        ops.integrator('LoadControl', 0.0)
        ok = ops.analyze(1)
        record()
        
        # Define integrator for main load
        ops.integrator('LoadControl', target_load_factor/steps)         

        # Run analysis
        for i in range(steps):
            if print_ops_status:
                print(f'Running Load Controlled Analysis Step {i}')

            ok = ops.analyze(1)

            if ok != 0:
                print(f'Load controlled analysis failed in step {i}')
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
                if results.absolute_maximum_strain[-1] > steel_strain_limit:
                    results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                    break
                    
            # Check for maximum PMM interaction value    
            if P_M_M_interaction_limit is not None:
                if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                    results.exit_message = 'P_M_M interaction Limit Reached'
                    break

        if not hasattr(results, 'exit_message'):
            results.exit_message = 'Full Load Applied'

        find_limit_point()
        return results



    def run_displacement_controlled_analysis(self, target_disp, steps=1000, **kwargs):
        steel_strain_limit = kwargs.get('steel_strain_limit', 0.05)
        eigenvalue_limit = kwargs.get('eigenvalue_limit', 0)
        P_M_M_interaction_limit=kwargs.get('P_M_M_interaction_limit', None)
        print_ops_status = kwargs.get('print_ops_status', True)

        # Define control node
        control_node = self.bottom_column_top_node_tag  # Bottom node of top column is the control node for displacement control
        print(f'Control node for displacement control: {control_node}')
        control_dof = 1   # Horizontal displacement is the control DOF

        # Initialize analysis results
        results = AnalysisResults(initialize_empty_lists = ['load_ratio',
                                                            'control_node_displacement',
                                                            'lowest_eigenvalue',
                                                            'vertical_reaction',
                                                            'lateral_reaction',
                                                            'absolute_maximum_strain',
                                                            'max_P_M_M_interaction',
                                                            ])

        # Define function to find limit point
        def find_limit_point():
            if 'Analysis Failed' in results.exit_message:
                ind, x = find_limit_point_in_list(results.load_ratio, max(results.load_ratio)) # @todo - should there be a limit point defined here?
            elif 'Eigenvalue Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.lowest_eigenvalue, eigenvalue_limit)
            elif 'Extreme Steel Fiber Strain Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.absolute_maximum_strain, steel_strain_limit)
            elif 'P_M_M interaction Limit Reached' in results.exit_message:
                ind, x = find_limit_point_in_list(results.max_P_M_M_interaction, P_M_M_interaction_limit)           
            else:
                raise Exception('Unknown limit point')
            
            results.maximum_load_ratio_at_limit_point = interpolate_list(results.load_ratio, ind, x)
            print(' Max Load Ratio',results.maximum_load_ratio_at_limit_point)

        # Define recorder
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
            results.absolute_maximum_strain.append(return_max_of_fiber_strain_in_all_elements(self.element_dict))
            max_PMM, max_ele_tag,_ = return_strength_ratio(self.element_dict)
            results.max_P_M_M_interaction.append(max_PMM)

        # Define analysis options
        ops.constraints('Transformation')
        ops.numberer('RCM')
        ops.system('UmfPack')
        ops.test('NormUnbalance', 1e-3, 10, 1)
        ops.algorithm('Newton')   
        ops.analysis('Static')
        
        # Run one step with no load
        ops.integrator('LoadControl', 0.0)
        ok = ops.analyze(1)
        record()

        # Define integrator for main load
        ops.integrator('DisplacementControl', control_node, control_dof, target_disp / steps)     

        # Run Displacement Control Analysis
        for i in range(steps):
            if print_ops_status:
                print(f'Running Displacement Controlled Analysis Step {i}')
            
            ok = ops.analyze(1)

            
            if ok != 0:
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
                print(results.absolute_maximum_strain)
                # input("Check the absolute maximum strain values. Press Enter to continue...")
                if results.absolute_maximum_strain[-1] > steel_strain_limit:
                    results.exit_message = 'Extreme Steel Fiber Strain Limit Reached'
                    break

            # Check for maximum PMM interaction value    
            if P_M_M_interaction_limit is not None:
                if results.max_P_M_M_interaction[-1] > P_M_M_interaction_limit:
                    results.exit_message = 'P_M_M interaction Limit Reached'
                    break

        if not hasattr(results, 'exit_message'):
            results.exit_message = 'Full Deformation Applied'
                              
        find_limit_point()
        return results

if __name__ == "__main__":
    kip = 1
    inch = 1
    ft = 12*inch
    ksi = kip/inch**2
    
    Stepped_Column = Stepped_Column(bottom_column_section_name='W14X132',
                                    height_of_bottom_column=8.0*ft,
                                    load_on_bottom_column=200*kip,
                                    top_column_section_name='W14X90',
                                    height_of_top_column=8.0*ft,
                                    number_of_elements=4,
                                    load_on_top_column=200*kip,
                                    lateral_load=20*kip,
                                    offset_1=0.3,
                                    offset_3=0.25,
                                    Fy=36*ksi,
                                    E=29000*ksi,
                                    Elastic_analysis=False,
                                    Second_order_effects=True,
                                    Residual_Stress=True,
                                    Geometric_Imperfection=True,
                                    geometric_imperfection_ratio=1/500,
                                    plot_model=False)

    Stepped_Column.build_stepped_column()
    #results= Stepped_Column.run_load_controlled_analysis(target_load_factor=3.0, steps=100)
    results = Stepped_Column.run_displacement_controlled_analysis(target_disp=1, steps=1000)
    print(results.exit_message)

    plot_deformed_2d(scale_factor=100,axis_equal=True)

    fig, ax = plt.subplots()
    plt.plot(results.load_ratio,results.control_node_displacement, marker = 'o', markersize=5)
    plt.xlabel('Load Ratio λ')
    plt.ylabel('Displacement at Control Node')

    fig, ax = plt.subplots()
    plt.plot(results.load_ratio,results.lowest_eigenvalue, marker = 'o', markersize=5)
    plt.xlabel('Load Ratio λ')
    plt.ylabel('Lowest Eigenvalue')
    
    fig, ax = plt.subplots()
    plt.plot(results.load_ratio,results.vertical_reaction, marker = 'o', markersize=5)
    plt.xlabel('Load Ratio λ')
    plt.ylabel('Vertical Reaction')
    
    fig, ax = plt.subplots()
    plt.plot(results.load_ratio,results.lateral_reaction, marker = 'o', markersize=5)
    plt.xlabel('Load Ratio λ')
    plt.ylabel('Lateral Reaction')
    
    fig, ax = plt.subplots()
    plt.plot(results.load_ratio,results.absolute_maximum_strain, marker = 'o', markersize=5)
    plt.xlabel('Load Ratio λ')
    plt.ylabel('Absolute Maximum Strain')
    
    fig, ax = plt.subplots()
    plt.plot(results.load_ratio,results.max_P_M_M_interaction, marker = 'o', markersize=5)
    plt.xlabel('Load Ratio λ')
    plt.ylabel('Max P-M-M Interaction')

    plt.show()