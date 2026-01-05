from Units import *


class convert_dict_items_to_class_attributes:
    def __init__(self,config):
        for k, v in config.items():
            setattr(self, k, v)    

Analysis_Info= {
    'GMNIA':{
                'Residual_Stress':True,
                'Elastic_analysis':False,
                'Second_order_effects':True,
                'stiffness_reduction':0.9,
                'strength_reduction':0.9,
                'Geometric_Imperfection':True,
                'Notional_load':False
                },
    'GMNA':{
                'Residual_Stress':True,
                'Elastic_analysis':False,
                'Second_order_effects':True,
                'stiffness_reduction':0.9,
                'strength_reduction':0.9,
                'Geometric_Imperfection':False,
                'Notional_load':False
                },
    'GNIA':{
                'Residual_Stress':False,
                'Elastic_analysis':True,
                'Second_order_effects':True,
                'stiffness_reduction':0.8,
                'strength_reduction':1,
                'Geometric_Imperfection':True,
                'Notional_load':False
                },
    'GNA_Notional_Loads':{
                'Residual_Stress':False,
                'Elastic_analysis':True,
                'Second_order_effects':True,
                'stiffness_reduction':0.8,
                'strength_reduction':1,
                'Geometric_Imperfection':False,
                'Notional_load':True
                },
    'GNA':{
                'Residual_Stress':False,
                'Elastic_analysis':True,
                'Second_order_effects':True,
                'stiffness_reduction':0.8,
                'strength_reduction':1,
                'Geometric_Imperfection':False,
                'Notional_load':False
                }
            }

Frame_Info={
        '0': {
        'Frame_id':'Test_Frame',
        'bay_width': [20 * ft, 48 * ft],
        'story_height': [20 * ft, 15 * ft],
        'column_no_of_ele': 1,
        'beam_no_of_ele': 1,
        'beam_section':
          {
            'common_and_exceptions': 
            {
                'common': 'W27X84',
                '(2,1)': 'W36X170',
                '(1,2)': 'W21X44',
                '(2,2)': 'W27X102',
             }
         },
        'column_section':
         {
            'common_and_exceptions': 
            {
                'common': ('W8X15', 'x'),
                '(2,1)': ('W14X132', 'y'),
                '(3,1)': ('W14X120', 'x'),
                '(1,2)': ('W8X13', 'y'),
                '(2,2)': ('W14X120', 'x'),
                '(3,2)': ('W14X109', 'y')
            }
        },
        'support': 'ppp',
        'load_comb_multipliers': [1.2, 1, 0.5, 1],
        'D_floor_intensity': 3.623 * kip / ft,
        'D_roof_intensity': 2.785 * kip / ft,
        'L_floor_intensity': 3.623 * kip / ft,
        'L_roof_intensity': 2.785 * kip / ft,
        'Wind_load_floor': 0,
        'Wind_load_roof': 0,
        'Wall_load':0,
        'geometric_imperfection_ratio': 1 / 500
         },

        '100': {
        'Frame_id':'Test_Frame',
        'bay_width': [20 * ft, 48 * ft],
        'story_height': [20 * ft, 15 * ft],
        'column_no_of_ele': 1,
        'beam_no_of_ele': 2,
        'beam_section':
          {
             'same_for_storey': 
            {
                '1': 'W14X30',
                '2': 'W14X26',
             }
         },
        'column_section':
         {
            'same_for_storey': {
                '1': ('W16X31','x'),
                '2': ('W16X31','y')
            }
        },
        'support': 'ppp',
        'load_comb_multipliers': [1.2, 1, 0.5, 1],
        'D_floor_intensity': 3.623 * kip / ft,
        'D_roof_intensity': 2.785 * kip / ft,
        'L_floor_intensity': 3.623 * kip / ft,
        'L_roof_intensity': 2.785 * kip / ft,
        'Wind_load_floor': 0,
        'Wind_load_roof': 0,
        'Wall_load':0,
        'geometric_imperfection_ratio': 1 / 500
         },
########################################################################## 


    '9': {
        'Frame_id':'Ziemian_9',
        'bay_width': [20 * ft, 48 * ft],
        'story_height': [20 * ft, 15 * ft],
        'column_no_of_ele': 4,
        'beam_no_of_ele': 4,
        'beam_section':
          {
            'common_and_exceptions': 
            {
                'common': 'W27X84',
                '(2,1)': 'W36X170',
                '(1,2)': 'W21X44',
                '(2,2)': 'W27X102',
             }
         },
        'column_section':
         {
            'common_and_exceptions': 
            {
                'common': ('W8X15', 'x'),
                '(2,1)': ('W14X132', 'x'),
                '(3,1)': ('W14X120', 'x'),
                '(1,2)': ('W8X13', 'x'),
                '(2,2)': ('W14X120', 'x'),
                '(3,2)': ('W14X109', 'x'),
            }
        },
        'support': 'ppp',
        'load_comb_multipliers': [1.2, 1, 0.5, 1],
        'D_floor_intensity': 3.623 * kip / ft,
        'D_roof_intensity': 2.785 * kip / ft,
        'L_floor_intensity': 3.623 * kip / ft,
        'L_roof_intensity': 2.785 * kip / ft,
        'Wind_load_floor': 0,
        'Wind_load_roof': 0,
        'Wall_load':0,
        'geometric_imperfection_ratio': 1 / 500
         },

######################################################################
    '10': {
        'Frame_id':'Ziemian_10',
        'bay_width': [20 * ft, 48 * ft],
        'story_height': [20 * ft, 15 * ft],
        'column_no_of_ele': 4,
        'beam_no_of_ele': 4,
        'beam_section':
          {
            'common_and_exceptions': 
            {
                'common': 'W30X99',
                '(2,1)': 'W36X160',
                '(1,2)': 'W24X55',
                '(2,2)': 'W30X99',
             }
         },
        'column_section':
         {
            'common_and_exceptions': 
            {
                'common': ('W14X90', 'y'),
                '(1,2)': ('W6X9', 'y'),
                '(2,2)': ('W8X40', 'y'),
                '(3,2)': ('W8X24', 'y'),
            }
        },
        'support': 'ppp',
        'load_comb_multipliers': [1.2, 1, 0.5, 1],
        'D_floor_intensity': 2.829 * kip / ft,
        'D_roof_intensity': 2.175 * kip / ft,
        'L_floor_intensity': 2.829 * kip / ft,
        'L_roof_intensity': 2.175 * kip / ft,
        'Wind_load_floor': 0,
        'Wind_load_roof': 0,
        'Wall_load':0,
        'geometric_imperfection_ratio': 1 / 500
         },

   ##########################################################################     
        
    '13': {
        'Frame_id':'Ziemian_13',
        'bay_width': [20 * ft, 20 * ft,20*ft],
        'story_height': [9*ft+6*inch]*10,
        'column_no_of_ele': 4,
        'beam_no_of_ele': 4,

        'beam_section' : {
            'same_for_storey': {
                '1': 'W16X31',
                '2': 'W16X31',
                '3': 'W14X30',
                '4': 'W14X30',
                '5': 'W14X26',
                '6': 'W14X26',
                '7': 'W14X22',
                '8': 'W14X22',
                '9': 'W14X22',
                '10': 'W12X22'
            }
        },
        'column_section' : {
            'common_and_exceptions': {
                '(1,1)': ('W8X67', 'x'),
                '(2,1)': ('W10X68', 'x'),
                '(3,1)': ('W10X68', 'x'),
                '(4,1)': ('W8X67', 'x'),

                '(1,2)': ('W8X67', 'x'),
                '(2,2)': ('W10X68', 'x'),
                '(3,2)': ('W10X68', 'x'),
                '(4,2)': ('W8X67', 'x'),

                '(1,3)': ('W8X58', 'x'),
                '(2,3)': ('W8X58', 'x'),
                '(3,3)': ('W8X58', 'x'),
                '(4,3)': ('W8X58', 'x'),

                '(1,4)': ('W8X58', 'x'),
                '(2,4)': ('W8X58', 'x'),
                '(3,4)': ('W8X58', 'x'),
                '(4,4)': ('W8X58', 'x'),

                '(1,5)': ('W8X48', 'x'),
                '(2,5)': ('W8X48', 'x'),
                '(3,5)': ('W8X48', 'x'),
                '(4,5)': ('W8X48', 'x'),

                '(1,6)': ('W8X48', 'x'),
                '(2,6)': ('W8X48', 'x'),
                '(3,6)': ('W8X48', 'x'),
                '(4,6)': ('W8X48', 'x'),

                '(1,7)': ('W8X35', 'x'),
                '(2,7)': ('W8X31', 'x'),
                '(3,7)': ('W8X31', 'x'),
                '(4,7)': ('W8X35', 'x'),

                '(1,8)': ('W8X35', 'x'),
                '(2,8)': ('W8X31', 'x'),
                '(3,8)': ('W8X31', 'x'),
                '(4,8)': ('W8X35', 'x'),

                '(1,9)': ('W8X24', 'x'),
                '(2,9)': ('W8X18', 'x'),
                '(3,9)': ('W8X18', 'x'),
                '(4,9)': ('W8X24', 'x'),

                '(1,10)': ('W8X24', 'x'),
                '(2,10)': ('W8X18', 'x'),
                '(3,10)': ('W8X18', 'x'),
                '(4,10)': ('W8X24', 'x')
            }
        },
        'support': 'ffff',
        'load_comb_multipliers': [1.2, 1, 0.5, 1],
        'D_floor_intensity': 1.069 * kip / ft,
        'D_roof_intensity': 1.034 * kip / ft,
        'L_floor_intensity': 0.778 * kip / ft,
        'L_roof_intensity': 0.778 * kip / ft,
        'Wind_load_floor': 6.379* kip,
        'Wind_load_roof': 3.189* kip,
        'Wall_load':9.234*kip,
        'geometric_imperfection_ratio': 1 / 500
         },

############## Ziemian Dissertation ###############################

        'UP36H': {
        'Frame_id':'Ziemian_UP36H',
        'bay_width': [20 * ft, 48 * ft],
        'story_height': [20 * ft, 15 * ft],
        'column_no_of_ele': 2,   ## change to 2
        'beam_no_of_ele': 4,     ## change to 4
        'beam_section':
          {
            'common_and_exceptions': 
            {
                'common': 'W27X84',
                '(2,1)': 'W36X170',
                '(1,2)': 'W21X44',
                '(2,2)': 'W27X102',
             }
         },
        'column_section':
         {
            'common_and_exceptions': 
            {
                'common': ('W8X15', 'x'),
                '(2,1)': ('W14X132', 'x'),
                '(3,1)': ('W14X120', 'x'),
                '(1,2)': ('W8X13', 'x'),
                '(2,2)': ('W14X120', 'x'),
                '(3,2)': ('W14X109', 'x'),
            }
        },
        'support': 'ppp',
        'load_comb_multipliers': [1.4, 0, 0, 0],
        'D_floor_intensity': 7.5 * kip / ft,
        'D_roof_intensity': 3.5 * kip / ft,
        'L_floor_intensity': 0 * kip / ft,
        'L_roof_intensity': 0 * kip / ft,
        'Wind_load_floor': 6.56,
        'Wind_load_roof': 2.81,
        'Wall_load':0,
        'geometric_imperfection_ratio': 1 / 500
         },


        'SP36H': {
        'Frame_id':'Ziemian_SP36H',
        'bay_width': [34 * ft, 34 * ft],
        'story_height': [20 * ft, 15 * ft],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section':
          {
            'common_and_exceptions': 
            {
                'common': 'W33X118',
                '(1,2)': 'W24X62',
                '(2,2)': 'W24X62',
             }
         },
        'column_section':
         {
            'common_and_exceptions': 
            {
                'common': ('W14X53', 'x'),
                '(2,1)': ('W14X90', 'x'),
                '(1,2)': ('W14X53', 'x'),
                '(2,2)': ('W14X22', 'x'),
                '(3,2)': ('W14X53', 'x'),
            }
        },
        'support': 'ppp',
        'load_comb_multipliers': [1.2, 0, 0, 1],
        'D_floor_intensity': 7.5 * kip / ft,
        'D_roof_intensity': 3.5 * kip / ft,
        'L_floor_intensity': 0 * kip / ft,
        'L_roof_intensity': 0 * kip / ft,
        'Wind_load_floor': 6.56*kip,
        'Wind_load_roof': 2.81*kip,
        'Wall_load':0,
        'geometric_imperfection_ratio': 1 / 500,
        'Leaning_column':True,
        'Leaning_column_offset': 2,
        'Leaning_column_floor_load':2,
        'Leaning_column_roof_load':1
         },



        'SF36H': {
        'Frame_id':'Ziemian_SF36H',
        'bay_width': [34 * ft, 34 * ft],
        'story_height': [20 * ft, 15 * ft],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section':
          {
            'common_and_exceptions': 
            {
                'common': 'W30X108',
                '(1,2)': 'W24X62',
                '(2,2)': 'W24X62',
             }
         },
        'column_section':
         {
            'common_and_exceptions': 
            {
                'common': ('W14X48', 'x'),
                '(2,1)': ('W14X74', 'x'),
                '(1,2)': ('W14X43', 'x'),
                '(2,2)': ('W14X43', 'x'),
                '(3,2)': ('W14X43', 'x'),
            }
        },
        'support': 'fff',
        'load_comb_multipliers': [1.2, 0, 0, 1],
        'D_floor_intensity': 7.5 * kip / ft,
        'D_roof_intensity': 3.5 * kip / ft,
        'L_floor_intensity': 0 * kip / ft,
        'L_roof_intensity': 0 * kip / ft,
        'Wind_load_floor': 6.56*kip,
        'Wind_load_roof': 2.81*kip,
        'Wall_load':0,
        'geometric_imperfection_ratio': 1 / 500
         },


        'SP36L': {
        'Frame_id':'Ziemian_SP36L',
        'bay_width': [34 * ft, 34 * ft],
        'story_height': [20 * ft, 15 * ft],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section':
          {
            'common_and_exceptions': 
            {
                'common': 'W24X55',
                '(1,2)': 'W16X31',
                '(2,2)': 'W16X31',
             }
         },
        'column_section':
         {
            'common_and_exceptions': 
            {
                'common': ('W14X61', 'x'),
                '(2,1)': ('W14X132', 'x'),
                '(1,2)': ('W12X14', 'x'),
                '(2,2)': ('W14X22', 'x'),
                '(3,2)': ('W12X14', 'x'),
            }
        },
        'support': 'ppp',
        'load_comb_multipliers': [1.4, 0, 0, 0],
        'D_floor_intensity': 2.25 * kip / ft,
        'D_roof_intensity': 1.125 * kip / ft,
        'L_floor_intensity': 0 * kip / ft,
        'L_roof_intensity': 0 * kip / ft,
        'Wind_load_floor': 6.56*kip,
        'Wind_load_roof': 2.81*kip,
        'Wall_load':0,
        'geometric_imperfection_ratio': 1 / 500
         },



        'UP36L': {
        'Frame_id':'Ziemian_UP36L',
        'bay_width': [20 * ft, 48 * ft],
        'story_height': [20 * ft, 15 * ft],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section':
          {
            'common_and_exceptions': 
            {
                'common': 'W24X62',
                '(2,1)': 'W24X68',
                '(1,2)': 'W12X19',
                '(2,2)': 'W21X44',
             }
         },
        'column_section':
         {
            'common_and_exceptions': 
            {
                'common': ('W14X61', 'x'),
                '(2,1)': ('W14X90', 'x'),
                '(3,1)': ('W14X48', 'x'),
                '(1,2)': ('W12X16', 'x'),
                '(2,2)': ('W14X38', 'x'),
                '(3,2)': ('W14X43', 'x'),
            }
        },
        'support': 'ppp',
        'load_comb_multipliers': [1.4, 0, 0, 0],
        'D_floor_intensity':2.25 * kip / ft,
        'D_roof_intensity': 1.125 * kip / ft,
        'L_floor_intensity': 0 * kip / ft,
        'L_roof_intensity': 0 * kip / ft,
        'Wind_load_floor': 6.56,
        'Wind_load_roof': 2.81,
        'Wall_load':0,
        'geometric_imperfection_ratio': 1 / 500
         }


     }
