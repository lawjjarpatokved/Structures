from helpers import ft,kip,inch,KN,m

Frame_Info={

###########################################################################

        'W14X132_x_1X12': {
            'Frame_id': 'W14X132_x_1X12',
            'bay_width': [],
            'story_height': [12 * ft],
            'column_no_of_ele': 2,
            'beam_no_of_ele': 4,
            'beam_section': {
                'common_and_exceptions': {
                    'common': 'W27X84',
                    '(2,1)': 'W36X170',
                    '(1,2)': 'W21X44',
                    '(2,2)': 'W27X102',
                }
            },
            'column_section': {
                'common_and_exceptions': {
                    'common': ('W14X132', 'x'),
                }
            },
            'support': 'f',
            'load_comb_multipliers': [1, 0, 0, 1],
            'D_floor_intensity': 7.5 * 10 * KN,
            'D_roof_intensity': 133.5 * KN,
            'L_floor_intensity': 0 * kip / ft,
            'L_roof_intensity': 0 * kip / ft,
            'Wind_load_floor': 6.56 * 3 * KN,
            'Wind_load_roof': 4.45 * 3 * KN,
            'Wall_load': 0,
            'geometric_imperfection_ratio': 1 / 500,
            'Leaning_column': False,
            'Leaning_column_offset': 2,
            'Leaning_column_floor_load': 2,
            'Leaning_column_roof_load': 1,
        },
############################################################################

        'W14X132_x_1X15': {
            'Frame_id': 'W14X132_x_1X15',
            'bay_width': [],
            'story_height': [15 * ft],
            'column_no_of_ele': 2,
            'beam_no_of_ele': 4,
            'beam_section': {
                'common_and_exceptions': {
                    'common': 'W27X84',
                    '(2,1)': 'W36X170',
                    '(1,2)': 'W21X44',
                    '(2,2)': 'W27X102',
                }
            },
            'column_section': {
                'common_and_exceptions': {
                    'common': ('W14X132', 'x'),
                }
            },
            'support': 'f',
            'load_comb_multipliers': [1, 0, 0, 1],
            'D_floor_intensity': 7.5 * 10 * KN,
            'D_roof_intensity': 133.5 * KN,
            'L_floor_intensity': 0 * kip / ft,
            'L_roof_intensity': 0 * kip / ft,
            'Wind_load_floor': 6.56 * 3 * KN,
            'Wind_load_roof': 4.45 * 3 * KN,
            'Wall_load': 0,
            'geometric_imperfection_ratio': 1 / 500,
            'Leaning_column': False,
            'Leaning_column_offset': 2,
            'Leaning_column_floor_load': 2,
            'Leaning_column_roof_load': 1,
        },

############################################################################

        'W14X132_x_1X18': {
            'Frame_id': 'W14X132_x_1X18',
            'bay_width': [],
            'story_height': [18 * ft],
            'column_no_of_ele': 2,
            'beam_no_of_ele': 4,
            'beam_section': {
                'common_and_exceptions': {
                    'common': 'W27X84',
                    '(2,1)': 'W36X170',
                    '(1,2)': 'W21X44',
                    '(2,2)': 'W27X102',
                }
            },
            'column_section': {
                'common_and_exceptions': {
                    'common': ('W14X132', 'x'),
                }
            },
            'support': 'f',
            'load_comb_multipliers': [1, 0, 0, 1],
            'D_floor_intensity': 7.5 * 10 * KN,
            'D_roof_intensity': 133.5 * KN,
            'L_floor_intensity': 0 * kip / ft,
            'L_roof_intensity': 0 * kip / ft,
            'Wind_load_floor': 6.56 * 3 * KN,
            'Wind_load_roof': 4.45 * 3 * KN,
            'Wall_load': 0,
            'geometric_imperfection_ratio': 1 / 500,
            'Leaning_column': False,
            'Leaning_column_offset': 2,
            'Leaning_column_floor_load': 2,
            'Leaning_column_roof_load': 1,
        },
#############################################################################
        'W14X132_x_1X30': {
            'Frame_id': 'W14X132_x_1X30',
            'bay_width': [],
            'story_height': [30 * ft],
            'column_no_of_ele': 2,
            'beam_no_of_ele': 4,
            'beam_section': {
                'common_and_exceptions': {
                    'common': 'W27X84',
                    '(2,1)': 'W36X170',
                    '(1,2)': 'W21X44',
                    '(2,2)': 'W27X102',
                }
            },
            'column_section': {
                'common_and_exceptions': {
                    'common': ('W14X132', 'x'),
                }
            },
            'support': 'f',
            'load_comb_multipliers': [1, 0, 0, 1],
            'D_floor_intensity': 7.5 * 10 * KN,
            'D_roof_intensity': 133.5 * KN,
            'L_floor_intensity': 0 * kip / ft,
            'L_roof_intensity': 0 * kip / ft,
            'Wind_load_floor': 6.56 * 3 * KN,
            'Wind_load_roof': 4.45 * 3 * KN,
            'Wall_load': 0,
            'geometric_imperfection_ratio': 1 / 500,
            'Leaning_column': False,
            'Leaning_column_offset': 2,
            'Leaning_column_floor_load': 2,
            'Leaning_column_roof_load': 1,
        },
#############################################################################
        'W14X132_x_1X36': {
            'Frame_id': 'W14X132_x_1X36',
            'bay_width': [],
            'story_height': [36 * ft],
            'column_no_of_ele': 2,
            'beam_no_of_ele': 4,
            'beam_section': {
                'common_and_exceptions': {
                    'common': 'W27X84',
                    '(2,1)': 'W36X170',
                    '(1,2)': 'W21X44',
                    '(2,2)': 'W27X102',
                }
            },
            'column_section': {
                'common_and_exceptions': {
                    'common': ('W14X132', 'x'),
                }
            },
            'support': 'f',
            'load_comb_multipliers': [1, 0, 0, 1],
            'D_floor_intensity': 7.5 * 10 * KN,
            'D_roof_intensity': 133.5 * KN,
            'L_floor_intensity': 0 * kip / ft,
            'L_roof_intensity': 0 * kip / ft,
            'Wind_load_floor': 6.56 * 3 * KN,
            'Wind_load_roof': 4.45 * 3 * KN,
            'Wall_load': 0,
            'geometric_imperfection_ratio': 1 / 500,
            'Leaning_column': False,
            'Leaning_column_offset': 2,
            'Leaning_column_floor_load': 2,
            'Leaning_column_roof_load': 1,
        },
#############################################################################
        'W27X84_x_1X12': {
            'Frame_id': 'W27X84_x_1X12',
            'bay_width': [],
            'story_height': [12 * ft],
            'column_no_of_ele': 2,
            'beam_no_of_ele': 4,
            'beam_section': {
                'common_and_exceptions': {
                    'common': 'W27X84',
                    '(2,1)': 'W36X170',
                    '(1,2)': 'W21X44',
                    '(2,2)': 'W27X102',
                }
            },
            'column_section': {
                'common_and_exceptions': {
                    'common': ('W27X84', 'x'),
                }
            },
            'support': 'f',
            'load_comb_multipliers': [1, 0, 0, 1],
            'D_floor_intensity': 7.5 * 10 * KN,
            'D_roof_intensity': 133.5 * KN,
            'L_floor_intensity': 0 * kip / ft,
            'L_roof_intensity': 0 * kip / ft,
            'Wind_load_floor': 6.56 * 3 * KN,
            'Wind_load_roof': 4.45 * 3 * KN,
            'Wall_load': 0,
            'geometric_imperfection_ratio': 1 / 500,
            'Leaning_column': False,
            'Leaning_column_offset': 2,
            'Leaning_column_floor_load': 2,
            'Leaning_column_roof_load': 1,
        },
############################################################################
#############################################################################
        'W27X84_x_1X18': {
            'Frame_id': 'W27X84_x_1X18',
            'bay_width': [],
            'story_height': [18 * ft],
            'column_no_of_ele': 2,
            'beam_no_of_ele': 4,
            'beam_section': {
                'common_and_exceptions': {
                    'common': 'W27X84',
                    '(2,1)': 'W36X170',
                    '(1,2)': 'W21X44',
                    '(2,2)': 'W27X102',
                }
            },
            'column_section': {
                'common_and_exceptions': {
                    'common': ('W27X84', 'x'),
                }
            },
            'support': 'f',
            'load_comb_multipliers': [1, 0, 0, 1],
            'D_floor_intensity': 7.5 * 10 * KN,
            'D_roof_intensity': 133.5 * KN,
            'L_floor_intensity': 0 * kip / ft,
            'L_roof_intensity': 0 * kip / ft,
            'Wind_load_floor': 6.56 * 3 * KN,
            'Wind_load_roof': 4.45 * 3 * KN,
            'Wall_load': 0,
            'geometric_imperfection_ratio': 1 / 500,
            'Leaning_column': False,
            'Leaning_column_offset': 2,
            'Leaning_column_floor_load': 2,
            'Leaning_column_roof_load': 1,
        },
############################################################################
#############################################################################
        'W27X84_x_1X30': {
            'Frame_id': 'W27X84_x_1X30',
            'bay_width': [],
            'story_height': [30 * ft],
            'column_no_of_ele': 2,
            'beam_no_of_ele': 4,
            'beam_section': {
                'common_and_exceptions': {
                    'common': 'W27X84',
                    '(2,1)': 'W36X170',
                    '(1,2)': 'W21X44',
                    '(2,2)': 'W27X102',
                }
            },
            'column_section': {
                'common_and_exceptions': {
                    'common': ('W27X84', 'x'),
                }
            },
            'support': 'f',
            'load_comb_multipliers': [1, 0, 0, 1],
            'D_floor_intensity': 7.5 * 10 * KN,
            'D_roof_intensity': 133.5 * KN,
            'L_floor_intensity': 0 * kip / ft,
            'L_roof_intensity': 0 * kip / ft,
            'Wind_load_floor': 6.56 * 3 * KN,
            'Wind_load_roof': 4.45 * 3 * KN,
            'Wall_load': 0,
            'geometric_imperfection_ratio': 1 / 500,
            'Leaning_column': False,
            'Leaning_column_offset': 2,
            'Leaning_column_floor_load': 2,
            'Leaning_column_roof_load': 1,
        },
############################################################################
        'Trial_Col': {
        'Frame_id':'Trial_Col',
        'bay_width': [],
        'story_height': [28 * ft],
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
            }
        },
        'support': 'f',
        'load_comb_multipliers': [1, 0, 0, 1],
        'D_floor_intensity': 7.5*10  * KN,
        'D_roof_intensity': 133.5 *KN ,
        'L_floor_intensity': 0 * kip / ft,
        'L_roof_intensity': 0 * kip / ft,
        'Wind_load_floor': 6.56*3*KN,
        'Wind_load_roof': 4.45*3*KN,
        'Wall_load':0,
        'geometric_imperfection_ratio':1/500,
        'Leaning_column':False,
        'Leaning_column_offset': 2,
        'Leaning_column_floor_load':2,
        'Leaning_column_roof_load':1,
         },
###########################################################################
        'Trial_Col_2': {
        'Frame_id':'Trial_Col_2',
        'bay_width': [],
        'story_height': [40 * ft, 28 * ft],
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
                'common': ('W14X132', 'x'),
                '(2,1)': ('W14X132', 'x'),
                '(3,1)': ('W14X120', 'x'),
                '(1,2)': ('W8X13', 'x'),
                '(2,2)': ('W14X120', 'x'),
                '(3,2)': ('W14X109', 'x'),
            }
        },
        'support': 'f',
        'load_comb_multipliers': [1, 0, 0, 1],
        'D_floor_intensity': 7.5*10  * KN,
        'D_roof_intensity': 133.5 *KN ,
        'L_floor_intensity': 0 * kip / ft,
        'L_roof_intensity': 0 * kip / ft,
        'Wind_load_floor': 6.56*3*KN,
        'Wind_load_roof': 4.45*3*KN,
        'Wall_load':0,
        'geometric_imperfection_ratio':1/500,
        'Leaning_column':False,
        'Leaning_column_offset': 2,
        'Leaning_column_floor_load':2,
        'Leaning_column_roof_load':1
         },
###########################################################################
        'W8X31_x_1X3': {
        'Frame_id': 'W8X31_x_1X3',
        'story_height': [3],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W8X31', 'x')}
          },
        'bay_width': [],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section':
          {
              'common_and_exceptions': {'common': 'W27X84'}
          },
        'support': 'f',
        'load_comb_multipliers': [1, 0, 0, 1],
        'D_floor_intensity': 75.0,
        'D_roof_intensity': 133.5,
        'L_floor_intensity': 0.0,
        'L_roof_intensity': 0.0,
        'Wind_load_floor': 19.68,
        'Wind_load_roof': 13.350000000000001,
        'Wall_load': 0,
        'geometric_imperfection_ratio': 0.002,
        'Leaning_column': False,
        'Leaning_column_offset': 2,
        'Leaning_column_floor_load': 2,
        'Leaning_column_roof_load': 1
        },
###########################################################################
        'W8X31_y_1X3': {
        'Frame_id': 'W8X31_y_1X3',
        'story_height': [3],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W8X31', 'y')}
          },
        'bay_width': [],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section':
          {
              'common_and_exceptions': {'common': 'W27X84'}
          },
        'support': 'f',
        'load_comb_multipliers': [1, 0, 0, 1],
        'D_floor_intensity': 75.0,
        'D_roof_intensity': 133.5,
        'L_floor_intensity': 0.0,
        'L_roof_intensity': 0.0,
        'Wind_load_floor': 19.68,
        'Wind_load_roof': 13.350000000000001,
        'Wall_load': 0,
        'geometric_imperfection_ratio': 0.002,
        'Leaning_column': False,
        'Leaning_column_offset': 2,
        'Leaning_column_floor_load': 2,
        'Leaning_column_roof_load': 1
        },
###########################################################################
        'W8X31_x_2X3': {
        'Frame_id': 'W8X31_x_2X3',
        'story_height': [3, 3],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W8X31', 'x')}
          },
        'bay_width': [],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section':
          {
              'common_and_exceptions': {'common': 'W27X84'}
          },
        'support': 'f',
        'load_comb_multipliers': [1, 0, 0, 1],
        'D_floor_intensity': 75.0,
        'D_roof_intensity': 133.5,
        'L_floor_intensity': 0.0,
        'L_roof_intensity': 0.0,
        'Wind_load_floor': 19.68,
        'Wind_load_roof': 13.350000000000001,
        'Wall_load': 0,
        'geometric_imperfection_ratio': 0.002,
        'Leaning_column': False,
        'Leaning_column_offset': 2,
        'Leaning_column_floor_load': 2,
        'Leaning_column_roof_load': 1
        },
###########################################################################
        'W8X31_y_2X3': {
        'Frame_id': 'W8X31_y_2X3',
        'story_height': [0.9143999999999999, 0.9143999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W8X31', 'y')}
          },
        'bay_width': [],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section':
          {
              'common_and_exceptions': {'common': 'W27X84'}
          },
        'support': 'f',
        'load_comb_multipliers': [1, 0, 0, 1],
        'D_floor_intensity': 75.0,
        'D_roof_intensity': 133.5,
        'L_floor_intensity': 0.0,
        'L_roof_intensity': 0.0,
        'Wind_load_floor': 19.68,
        'Wind_load_roof': 13.350000000000001,
        'Wall_load': 0,
        'geometric_imperfection_ratio': 0.002,
        'Leaning_column': False,
        'Leaning_column_offset': 2,
        'Leaning_column_floor_load': 2,
        'Leaning_column_roof_load': 1
        },
###########################################################################
        'W8X31_x_1X4': {
        'Frame_id': 'W8X31_x_1X4',
        'story_height': [1.2191999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W8X31', 'x')}
          },
        'bay_width': [],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section':
          {
              'common_and_exceptions': {'common': 'W27X84'}
          },
        'support': 'f',
        'load_comb_multipliers': [1, 0, 0, 1],
        'D_floor_intensity': 75.0,
        'D_roof_intensity': 133.5,
        'L_floor_intensity': 0.0,
        'L_roof_intensity': 0.0,
        'Wind_load_floor': 19.68,
        'Wind_load_roof': 13.350000000000001,
        'Wall_load': 0,
        'geometric_imperfection_ratio': 0.002,
        'Leaning_column': False,
        'Leaning_column_offset': 2,
        'Leaning_column_floor_load': 2,
        'Leaning_column_roof_load': 1
        }
}
