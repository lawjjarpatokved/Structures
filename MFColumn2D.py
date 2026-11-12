from Structures_2D import Structures_2D



class MFColumn_2D(Structures_2D):

    print_ops_status=True

    def __init__(self,storey_height,
                no_of_elements_column,
                 column_section,load_combination_multipliers,Frame_id,Material_obj,Steel_Grade,
                 **kwargs):
        width_of_bay=[]
        no_of_elements_beam= 4
        beam_section= {
            "common_and_exceptions": {
                "common": "W27X84"
            }
        }
        kwargs["column_only_model"] = True 
        super().__init__(width_of_bay,storey_height,
                no_of_elements_column, no_of_elements_beam,
                 beam_section,column_section,load_combination_multipliers,Frame_id,Material_obj,Steel_Grade,
                 **kwargs)

    @classmethod
    def from_dict(cls, details_dict: dict):
        return cls(**details_dict)

    def axis_i_floor_load(self, i=1):
        if i == 1:
            return [(-self.D_multiplier*self.D_floor_intensity)+(-self.L_multiplier*self.L_floor_intensity)]
        else:
            raise Exception("axis_i_floor_load is only defined for i=1 in MFColumn_2D class as there is only one axis in column models.")

    def axis_i_roof_load(self, i=1):
        if i == 1:
            return [(-self.D_multiplier*self.D_roof_intensity)+(-self.L_r_multiplier*self.L_roof_intensity)]
        else:
            raise Exception("axis_i_roof_load is only defined for i=1 in MFColumn_2D class as there is only one axis in column models.")
        
    def floor_notional_load(self):
        return (((self.D_multiplier*self.D_floor_intensity)+(self.L_multiplier*self.L_floor_intensity)))*0.002

    def roof_notional_load(self):
        return (((self.D_multiplier*self.D_roof_intensity)+(self.L_r_multiplier*self.L_roof_intensity)))*0.002
    
    def get_control_node_and_dof(self,control_dir='L'):
        if control_dir=='L':
            control_node = self.axis_i_roof_nodes(1)[0]   
            control_dof = 1  # horizontal displacement
        # Control node for vertical deflection of beam node
        else:
            control_node = self.axis_i_roof_nodes(1)[0]  
            control_dof = 2  # vertical displacement    

        return control_node,control_dof    
