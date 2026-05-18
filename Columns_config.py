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
        },
###########################################################################
        'W8X31_x_1X0.9143999999999999': {
        'Frame_id': 'W8X31_x_1X0.9143999999999999',
        'story_height': [0.9143999999999999],
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
        'W8X31_x_1X1.8287999999999998': {
        'Frame_id': 'W8X31_x_1X1.8287999999999998',
        'story_height': [1.8287999999999998],
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
        'W8X31_x_1X2.1335999999999995': {
        'Frame_id': 'W8X31_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
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
        'W8X31_x_1X2.2859999999999996': {
        'Frame_id': 'W8X31_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
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
        'W8X31_x_1X2.4383999999999997': {
        'Frame_id': 'W8X31_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
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
        'W8X31_x_1X2.5907999999999998': {
        'Frame_id': 'W8X31_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
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
        'W8X31_x_1X2.7432': {
        'Frame_id': 'W8X31_x_1X2.7432',
        'story_height': [2.7432],
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
        'W8X31_x_1X2.8955999999999995': {
        'Frame_id': 'W8X31_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
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
        'W8X31_x_1X3.0479999999999996': {
        'Frame_id': 'W8X31_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
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
        'W8X31_x_1X3.2003999999999997': {
        'Frame_id': 'W8X31_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
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
        'W8X31_x_1X3.3527999999999993': {
        'Frame_id': 'W8X31_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
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
        'W8X31_x_1X3.5051999999999994': {
        'Frame_id': 'W8X31_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
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
        'W8X31_x_1X3.6575999999999995': {
        'Frame_id': 'W8X31_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
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
        'W8X31_x_1X3.8099999999999996': {
        'Frame_id': 'W8X31_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
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
        'W8X31_x_1X3.9623999999999997': {
        'Frame_id': 'W8X31_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
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
        'W8X31_x_1X4.1148': {
        'Frame_id': 'W8X31_x_1X4.1148',
        'story_height': [4.1148],
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
        'W8X31_x_1X4.267199999999999': {
        'Frame_id': 'W8X31_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
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
        'W8X31_x_1X4.419599999999999': {
        'Frame_id': 'W8X31_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
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
        'W8X31_x_1X4.571999999999999': {
        'Frame_id': 'W8X31_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
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
        'W8X31_x_1X4.724399999999999': {
        'Frame_id': 'W8X31_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
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
        'W8X31_x_1X4.876799999999999': {
        'Frame_id': 'W8X31_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
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
        'W8X31_x_1X5.0291999999999994': {
        'Frame_id': 'W8X31_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
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
        'W8X31_x_1X5.1815999999999995': {
        'Frame_id': 'W8X31_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
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
        'W8X31_x_1X5.334': {
        'Frame_id': 'W8X31_x_1X5.334',
        'story_height': [5.334],
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
        'W8X31_x_1X5.4864': {
        'Frame_id': 'W8X31_x_1X5.4864',
        'story_height': [5.4864],
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
        'W8X31_x_1X5.638799999999999': {
        'Frame_id': 'W8X31_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
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
        'W8X31_x_1X5.791199999999999': {
        'Frame_id': 'W8X31_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
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
        'W8X31_x_1X5.943599999999999': {
        'Frame_id': 'W8X31_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
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
        'W8X31_x_1X6.095999999999999': {
        'Frame_id': 'W8X31_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
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
        'W14X43_x_1X2.1335999999999995': {
        'Frame_id': 'W14X43_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X2.2859999999999996': {
        'Frame_id': 'W14X43_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X2.4383999999999997': {
        'Frame_id': 'W14X43_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X2.5907999999999998': {
        'Frame_id': 'W14X43_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X2.7432': {
        'Frame_id': 'W14X43_x_1X2.7432',
        'story_height': [2.7432],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X2.8955999999999995': {
        'Frame_id': 'W14X43_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X3.0479999999999996': {
        'Frame_id': 'W14X43_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X3.2003999999999997': {
        'Frame_id': 'W14X43_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X3.3527999999999993': {
        'Frame_id': 'W14X43_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X3.5051999999999994': {
        'Frame_id': 'W14X43_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X3.6575999999999995': {
        'Frame_id': 'W14X43_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X3.8099999999999996': {
        'Frame_id': 'W14X43_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X3.9623999999999997': {
        'Frame_id': 'W14X43_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X4.1148': {
        'Frame_id': 'W14X43_x_1X4.1148',
        'story_height': [4.1148],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X4.267199999999999': {
        'Frame_id': 'W14X43_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X4.419599999999999': {
        'Frame_id': 'W14X43_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X4.571999999999999': {
        'Frame_id': 'W14X43_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X4.724399999999999': {
        'Frame_id': 'W14X43_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X4.876799999999999': {
        'Frame_id': 'W14X43_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X5.0291999999999994': {
        'Frame_id': 'W14X43_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X5.1815999999999995': {
        'Frame_id': 'W14X43_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X5.334': {
        'Frame_id': 'W14X43_x_1X5.334',
        'story_height': [5.334],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X5.4864': {
        'Frame_id': 'W14X43_x_1X5.4864',
        'story_height': [5.4864],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X5.638799999999999': {
        'Frame_id': 'W14X43_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X5.791199999999999': {
        'Frame_id': 'W14X43_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X5.943599999999999': {
        'Frame_id': 'W14X43_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X43_x_1X6.095999999999999': {
        'Frame_id': 'W14X43_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X43', 'x')}
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
        'W14X120_x_1X2.1335999999999995': {
        'Frame_id': 'W14X120_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X2.2859999999999996': {
        'Frame_id': 'W14X120_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X2.4383999999999997': {
        'Frame_id': 'W14X120_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X2.5907999999999998': {
        'Frame_id': 'W14X120_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X2.7432': {
        'Frame_id': 'W14X120_x_1X2.7432',
        'story_height': [2.7432],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X2.8955999999999995': {
        'Frame_id': 'W14X120_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X3.0479999999999996': {
        'Frame_id': 'W14X120_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X3.2003999999999997': {
        'Frame_id': 'W14X120_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X3.3527999999999993': {
        'Frame_id': 'W14X120_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X3.5051999999999994': {
        'Frame_id': 'W14X120_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X3.6575999999999995': {
        'Frame_id': 'W14X120_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X3.8099999999999996': {
        'Frame_id': 'W14X120_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X3.9623999999999997': {
        'Frame_id': 'W14X120_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X4.1148': {
        'Frame_id': 'W14X120_x_1X4.1148',
        'story_height': [4.1148],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X4.267199999999999': {
        'Frame_id': 'W14X120_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X4.419599999999999': {
        'Frame_id': 'W14X120_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X4.571999999999999': {
        'Frame_id': 'W14X120_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X4.724399999999999': {
        'Frame_id': 'W14X120_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X4.876799999999999': {
        'Frame_id': 'W14X120_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X5.0291999999999994': {
        'Frame_id': 'W14X120_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X5.1815999999999995': {
        'Frame_id': 'W14X120_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X5.334': {
        'Frame_id': 'W14X120_x_1X5.334',
        'story_height': [5.334],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X5.4864': {
        'Frame_id': 'W14X120_x_1X5.4864',
        'story_height': [5.4864],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X5.638799999999999': {
        'Frame_id': 'W14X120_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X5.791199999999999': {
        'Frame_id': 'W14X120_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X5.943599999999999': {
        'Frame_id': 'W14X120_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X120_x_1X6.095999999999999': {
        'Frame_id': 'W14X120_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X120', 'x')}
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
        'W14X311_x_1X2.1335999999999995': {
        'Frame_id': 'W14X311_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X2.2859999999999996': {
        'Frame_id': 'W14X311_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X2.4383999999999997': {
        'Frame_id': 'W14X311_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X2.5907999999999998': {
        'Frame_id': 'W14X311_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X2.7432': {
        'Frame_id': 'W14X311_x_1X2.7432',
        'story_height': [2.7432],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X2.8955999999999995': {
        'Frame_id': 'W14X311_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X3.0479999999999996': {
        'Frame_id': 'W14X311_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X3.2003999999999997': {
        'Frame_id': 'W14X311_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X3.3527999999999993': {
        'Frame_id': 'W14X311_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X3.5051999999999994': {
        'Frame_id': 'W14X311_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X3.6575999999999995': {
        'Frame_id': 'W14X311_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X3.8099999999999996': {
        'Frame_id': 'W14X311_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X3.9623999999999997': {
        'Frame_id': 'W14X311_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X4.1148': {
        'Frame_id': 'W14X311_x_1X4.1148',
        'story_height': [4.1148],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X4.267199999999999': {
        'Frame_id': 'W14X311_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X4.419599999999999': {
        'Frame_id': 'W14X311_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X4.571999999999999': {
        'Frame_id': 'W14X311_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X4.724399999999999': {
        'Frame_id': 'W14X311_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X4.876799999999999': {
        'Frame_id': 'W14X311_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X5.0291999999999994': {
        'Frame_id': 'W14X311_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X5.1815999999999995': {
        'Frame_id': 'W14X311_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X5.334': {
        'Frame_id': 'W14X311_x_1X5.334',
        'story_height': [5.334],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X5.4864': {
        'Frame_id': 'W14X311_x_1X5.4864',
        'story_height': [5.4864],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X5.638799999999999': {
        'Frame_id': 'W14X311_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X5.791199999999999': {
        'Frame_id': 'W14X311_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X5.943599999999999': {
        'Frame_id': 'W14X311_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W14X311_x_1X6.095999999999999': {
        'Frame_id': 'W14X311_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X311', 'x')}
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
        'W21X62_x_1X2.1335999999999995': {
        'Frame_id': 'W21X62_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X2.2859999999999996': {
        'Frame_id': 'W21X62_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X2.4383999999999997': {
        'Frame_id': 'W21X62_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X2.5907999999999998': {
        'Frame_id': 'W21X62_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X2.7432': {
        'Frame_id': 'W21X62_x_1X2.7432',
        'story_height': [2.7432],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X2.8955999999999995': {
        'Frame_id': 'W21X62_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X3.0479999999999996': {
        'Frame_id': 'W21X62_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X3.2003999999999997': {
        'Frame_id': 'W21X62_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X3.3527999999999993': {
        'Frame_id': 'W21X62_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X3.5051999999999994': {
        'Frame_id': 'W21X62_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X3.6575999999999995': {
        'Frame_id': 'W21X62_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X3.8099999999999996': {
        'Frame_id': 'W21X62_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X3.9623999999999997': {
        'Frame_id': 'W21X62_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X4.1148': {
        'Frame_id': 'W21X62_x_1X4.1148',
        'story_height': [4.1148],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X4.267199999999999': {
        'Frame_id': 'W21X62_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X4.419599999999999': {
        'Frame_id': 'W21X62_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X4.571999999999999': {
        'Frame_id': 'W21X62_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X4.724399999999999': {
        'Frame_id': 'W21X62_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X4.876799999999999': {
        'Frame_id': 'W21X62_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X5.0291999999999994': {
        'Frame_id': 'W21X62_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X5.1815999999999995': {
        'Frame_id': 'W21X62_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X5.334': {
        'Frame_id': 'W21X62_x_1X5.334',
        'story_height': [5.334],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X5.4864': {
        'Frame_id': 'W21X62_x_1X5.4864',
        'story_height': [5.4864],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X5.638799999999999': {
        'Frame_id': 'W21X62_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X5.791199999999999': {
        'Frame_id': 'W21X62_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X5.943599999999999': {
        'Frame_id': 'W21X62_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W21X62_x_1X6.095999999999999': {
        'Frame_id': 'W21X62_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W21X62', 'x')}
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
        'W40X392_x_1X2.1335999999999995': {
        'Frame_id': 'W40X392_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X2.2859999999999996': {
        'Frame_id': 'W40X392_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X2.4383999999999997': {
        'Frame_id': 'W40X392_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X2.5907999999999998': {
        'Frame_id': 'W40X392_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X2.7432': {
        'Frame_id': 'W40X392_x_1X2.7432',
        'story_height': [2.7432],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X2.8955999999999995': {
        'Frame_id': 'W40X392_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X3.0479999999999996': {
        'Frame_id': 'W40X392_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X3.2003999999999997': {
        'Frame_id': 'W40X392_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X3.3527999999999993': {
        'Frame_id': 'W40X392_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X3.5051999999999994': {
        'Frame_id': 'W40X392_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X3.6575999999999995': {
        'Frame_id': 'W40X392_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X3.8099999999999996': {
        'Frame_id': 'W40X392_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X3.9623999999999997': {
        'Frame_id': 'W40X392_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X4.1148': {
        'Frame_id': 'W40X392_x_1X4.1148',
        'story_height': [4.1148],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X4.267199999999999': {
        'Frame_id': 'W40X392_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X4.419599999999999': {
        'Frame_id': 'W40X392_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X4.571999999999999': {
        'Frame_id': 'W40X392_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X4.724399999999999': {
        'Frame_id': 'W40X392_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X4.876799999999999': {
        'Frame_id': 'W40X392_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X5.0291999999999994': {
        'Frame_id': 'W40X392_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X5.1815999999999995': {
        'Frame_id': 'W40X392_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X5.334': {
        'Frame_id': 'W40X392_x_1X5.334',
        'story_height': [5.334],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X5.4864': {
        'Frame_id': 'W40X392_x_1X5.4864',
        'story_height': [5.4864],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X5.638799999999999': {
        'Frame_id': 'W40X392_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X5.791199999999999': {
        'Frame_id': 'W40X392_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X5.943599999999999': {
        'Frame_id': 'W40X392_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W40X392_x_1X6.095999999999999': {
        'Frame_id': 'W40X392_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X392', 'x')}
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
        'W14X730_x_1X2.1335999999999995': {
        'Frame_id': 'W14X730_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X2.2859999999999996': {
        'Frame_id': 'W14X730_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X2.4383999999999997': {
        'Frame_id': 'W14X730_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X2.5907999999999998': {
        'Frame_id': 'W14X730_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X2.7432': {
        'Frame_id': 'W14X730_x_1X2.7432',
        'story_height': [2.7432],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X2.8955999999999995': {
        'Frame_id': 'W14X730_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X3.0479999999999996': {
        'Frame_id': 'W14X730_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X3.2003999999999997': {
        'Frame_id': 'W14X730_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X3.3527999999999993': {
        'Frame_id': 'W14X730_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X3.5051999999999994': {
        'Frame_id': 'W14X730_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X3.6575999999999995': {
        'Frame_id': 'W14X730_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X3.8099999999999996': {
        'Frame_id': 'W14X730_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X3.9623999999999997': {
        'Frame_id': 'W14X730_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X4.1148': {
        'Frame_id': 'W14X730_x_1X4.1148',
        'story_height': [4.1148],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X4.267199999999999': {
        'Frame_id': 'W14X730_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X4.419599999999999': {
        'Frame_id': 'W14X730_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X4.571999999999999': {
        'Frame_id': 'W14X730_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X4.724399999999999': {
        'Frame_id': 'W14X730_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X4.876799999999999': {
        'Frame_id': 'W14X730_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X5.0291999999999994': {
        'Frame_id': 'W14X730_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X5.1815999999999995': {
        'Frame_id': 'W14X730_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X5.334': {
        'Frame_id': 'W14X730_x_1X5.334',
        'story_height': [5.334],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X5.4864': {
        'Frame_id': 'W14X730_x_1X5.4864',
        'story_height': [5.4864],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X5.638799999999999': {
        'Frame_id': 'W14X730_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X5.791199999999999': {
        'Frame_id': 'W14X730_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X5.943599999999999': {
        'Frame_id': 'W14X730_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W14X730_x_1X6.095999999999999': {
        'Frame_id': 'W14X730_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W14X730', 'x')}
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
        'W18X311_x_1X2.1335999999999995': {
        'Frame_id': 'W18X311_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X2.2859999999999996': {
        'Frame_id': 'W18X311_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X2.4383999999999997': {
        'Frame_id': 'W18X311_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X2.5907999999999998': {
        'Frame_id': 'W18X311_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X2.7432': {
        'Frame_id': 'W18X311_x_1X2.7432',
        'story_height': [2.7432],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X2.8955999999999995': {
        'Frame_id': 'W18X311_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X3.0479999999999996': {
        'Frame_id': 'W18X311_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X3.2003999999999997': {
        'Frame_id': 'W18X311_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X3.3527999999999993': {
        'Frame_id': 'W18X311_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X3.5051999999999994': {
        'Frame_id': 'W18X311_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X3.6575999999999995': {
        'Frame_id': 'W18X311_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X3.8099999999999996': {
        'Frame_id': 'W18X311_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X3.9623999999999997': {
        'Frame_id': 'W18X311_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X4.1148': {
        'Frame_id': 'W18X311_x_1X4.1148',
        'story_height': [4.1148],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X4.267199999999999': {
        'Frame_id': 'W18X311_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X4.419599999999999': {
        'Frame_id': 'W18X311_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X4.571999999999999': {
        'Frame_id': 'W18X311_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X4.724399999999999': {
        'Frame_id': 'W18X311_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X4.876799999999999': {
        'Frame_id': 'W18X311_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X5.0291999999999994': {
        'Frame_id': 'W18X311_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X5.1815999999999995': {
        'Frame_id': 'W18X311_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X5.334': {
        'Frame_id': 'W18X311_x_1X5.334',
        'story_height': [5.334],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X5.4864': {
        'Frame_id': 'W18X311_x_1X5.4864',
        'story_height': [5.4864],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X5.638799999999999': {
        'Frame_id': 'W18X311_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X5.791199999999999': {
        'Frame_id': 'W18X311_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X5.943599999999999': {
        'Frame_id': 'W18X311_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W18X311_x_1X6.095999999999999': {
        'Frame_id': 'W18X311_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W18X311', 'x')}
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
        'W40X264_x_1X2.1335999999999995': {
        'Frame_id': 'W40X264_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X2.2859999999999996': {
        'Frame_id': 'W40X264_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X2.4383999999999997': {
        'Frame_id': 'W40X264_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X2.5907999999999998': {
        'Frame_id': 'W40X264_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X2.7432': {
        'Frame_id': 'W40X264_x_1X2.7432',
        'story_height': [2.7432],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X2.8955999999999995': {
        'Frame_id': 'W40X264_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X3.0479999999999996': {
        'Frame_id': 'W40X264_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X3.2003999999999997': {
        'Frame_id': 'W40X264_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X3.3527999999999993': {
        'Frame_id': 'W40X264_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X3.5051999999999994': {
        'Frame_id': 'W40X264_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X3.6575999999999995': {
        'Frame_id': 'W40X264_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X3.8099999999999996': {
        'Frame_id': 'W40X264_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X3.9623999999999997': {
        'Frame_id': 'W40X264_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X4.1148': {
        'Frame_id': 'W40X264_x_1X4.1148',
        'story_height': [4.1148],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X4.267199999999999': {
        'Frame_id': 'W40X264_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X4.419599999999999': {
        'Frame_id': 'W40X264_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X4.571999999999999': {
        'Frame_id': 'W40X264_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X4.724399999999999': {
        'Frame_id': 'W40X264_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X4.876799999999999': {
        'Frame_id': 'W40X264_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X5.0291999999999994': {
        'Frame_id': 'W40X264_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X5.1815999999999995': {
        'Frame_id': 'W40X264_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X5.334': {
        'Frame_id': 'W40X264_x_1X5.334',
        'story_height': [5.334],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X5.4864': {
        'Frame_id': 'W40X264_x_1X5.4864',
        'story_height': [5.4864],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X5.638799999999999': {
        'Frame_id': 'W40X264_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X5.791199999999999': {
        'Frame_id': 'W40X264_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X5.943599999999999': {
        'Frame_id': 'W40X264_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X264_x_1X6.095999999999999': {
        'Frame_id': 'W40X264_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X264', 'x')}
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
        'W40X593_x_1X2.1335999999999995': {
        'Frame_id': 'W40X593_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X2.2859999999999996': {
        'Frame_id': 'W40X593_x_1X2.2859999999999996',
        'story_height': [2.2859999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X2.4383999999999997': {
        'Frame_id': 'W40X593_x_1X2.4383999999999997',
        'story_height': [2.4383999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X2.5907999999999998': {
        'Frame_id': 'W40X593_x_1X2.5907999999999998',
        'story_height': [2.5907999999999998],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X2.7432': {
        'Frame_id': 'W40X593_x_1X2.7432',
        'story_height': [2.7432],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X2.8955999999999995': {
        'Frame_id': 'W40X593_x_1X2.8955999999999995',
        'story_height': [2.8955999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X3.0479999999999996': {
        'Frame_id': 'W40X593_x_1X3.0479999999999996',
        'story_height': [3.0479999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X3.2003999999999997': {
        'Frame_id': 'W40X593_x_1X3.2003999999999997',
        'story_height': [3.2003999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X3.3527999999999993': {
        'Frame_id': 'W40X593_x_1X3.3527999999999993',
        'story_height': [3.3527999999999993],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X3.5051999999999994': {
        'Frame_id': 'W40X593_x_1X3.5051999999999994',
        'story_height': [3.5051999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X3.6575999999999995': {
        'Frame_id': 'W40X593_x_1X3.6575999999999995',
        'story_height': [3.6575999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X3.8099999999999996': {
        'Frame_id': 'W40X593_x_1X3.8099999999999996',
        'story_height': [3.8099999999999996],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X3.9623999999999997': {
        'Frame_id': 'W40X593_x_1X3.9623999999999997',
        'story_height': [3.9623999999999997],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X4.1148': {
        'Frame_id': 'W40X593_x_1X4.1148',
        'story_height': [4.1148],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X4.267199999999999': {
        'Frame_id': 'W40X593_x_1X4.267199999999999',
        'story_height': [4.267199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X4.419599999999999': {
        'Frame_id': 'W40X593_x_1X4.419599999999999',
        'story_height': [4.419599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X4.571999999999999': {
        'Frame_id': 'W40X593_x_1X4.571999999999999',
        'story_height': [4.571999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X4.724399999999999': {
        'Frame_id': 'W40X593_x_1X4.724399999999999',
        'story_height': [4.724399999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X4.876799999999999': {
        'Frame_id': 'W40X593_x_1X4.876799999999999',
        'story_height': [4.876799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X5.0291999999999994': {
        'Frame_id': 'W40X593_x_1X5.0291999999999994',
        'story_height': [5.0291999999999994],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X5.1815999999999995': {
        'Frame_id': 'W40X593_x_1X5.1815999999999995',
        'story_height': [5.1815999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X5.334': {
        'Frame_id': 'W40X593_x_1X5.334',
        'story_height': [5.334],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X5.4864': {
        'Frame_id': 'W40X593_x_1X5.4864',
        'story_height': [5.4864],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X5.638799999999999': {
        'Frame_id': 'W40X593_x_1X5.638799999999999',
        'story_height': [5.638799999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X5.791199999999999': {
        'Frame_id': 'W40X593_x_1X5.791199999999999',
        'story_height': [5.791199999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X5.943599999999999': {
        'Frame_id': 'W40X593_x_1X5.943599999999999',
        'story_height': [5.943599999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        'W40X593_x_1X6.095999999999999': {
        'Frame_id': 'W40X593_x_1X6.095999999999999',
        'story_height': [6.095999999999999],
        'column_section':
          {
              'common_and_exceptions': {'common': ('W40X593', 'x')}
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
        '8X31_x_1X2.1335999999999995': {
        'Frame_id': '8X31_x_1X2.1335999999999995',
        'story_height': [2.1335999999999995],
        'column_section':
          {
              'common_and_exceptions': {'common': ('8X31', 'x')}
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
