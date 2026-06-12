m=1.0; mm=0.001*m; cm=0.01*m; km=1000.0*m; inch= 0.0254*m;ft=inch*12
KN=1.0; N=0.001*KN
Mpa= 10**3* KN/m**2; Gpa= 10**6* KN/m**2
sec=1
tonne=1*KN*(sec**2)/m; kg=0.001*tonne
lb = 4.4482216152605 * N   # 1 lbf = 4.4482216152605 N
kip = 1000.0 * lb
slug = lb * sec**2 / ft    # mass in US (F = m·a)
lbm = 0.45359237 * kg      # 1 lbm = 0.45359237 kg
psi = lb / inch**2
ksi = kip / inch**2

import libdenavit.section.database.aisc as section
import libdenavit.section.wide_flange as database
import os
import json
from Columns_config import Frame_Info
from Analysis_config import Analysis_Info
from Materials_config import Material_Info

class wf_Database:
    def __init__(self,Section_name,unit=inch):

        self.section=Section_name
        self.d=section.wide_flange_database[self.section]['d']*unit
        self.tw=section.wide_flange_database[self.section]['tw']*unit
        self.bf=section.wide_flange_database[self.section]['bf']*unit
        self.tf=section.wide_flange_database[self.section]['tf']*unit
        self.A=section.wide_flange_database[self.section]['A']*(unit**2)
        self.Ix=section.wide_flange_database[self.section]['Ix']*(unit**4)
        self.Iy=section.wide_flange_database[self.section]['Iy']*(unit**4)

class WF_Database:
        def __init__(self,Section_name,unit=inch):
        
            db = database.WideFlangeDB(Section_name)
            self.d   = db.d   * unit
            self.tw  = db.tw  * unit
            self.bf  = db.bf  * unit
            self.tf  = db.tf  * unit
            self.A   = db.A   * unit**2
            self.Ix  = db.Ix  * unit**4
            self.Zx  = db.Zx  * unit**3
            self.Sx  = db.Sx  * unit**3
            self.rx  = db.rx  * unit
            self.Iy  = db.Iy  * unit**4
            self.Zy  = db.Zy  * unit**3
            self.Sy  = db.Sy  * unit**3
            self.ry  = db.ry  * unit
            self.J   = db.J   * unit**4
            self.Cw  = db.Cw  * unit**6
            self.rts = db.rts * unit
            self.ho  = db.ho  * unit

class Steel_Material:
    def __init__(self,mat_tag,E,Fy):
        self.mat_tag=mat_tag
        self.E=E
        self.Fy=Fy

        
class convert_dict_items_to_class_attributes:
    def __init__(self,config):
        for k, v in config.items():
            setattr(self, k, v)    

def load_wind_dirn_data(json_wind_dirn_path):
    if not os.path.exists(json_wind_dirn_path):
        return {}

    try:
        with open(json_wind_dirn_path, "r") as f:
            content = f.read().strip()

            if content == "":
                return {}

            return json.loads(content)

    except json.JSONDecodeError:
        return {}

def save_wind_dirn_data(data,json_wind_dirn_path):
    with open(json_wind_dirn_path, "w") as f:
        json.dump(data, f, indent=4)

def ensure_frame_entry_exists(frame_key,data,json_wind_dirn_path):
    if frame_key not in data:
        data[frame_key]={
            "wind_load_dirn": None,
            "wind_load_dirn_source": "unknown"
        }

        save_wind_dirn_data(data,json_wind_dirn_path=json_wind_dirn_path)
    return data

def _format_frame_entry(code, config):
    """Format a frame configuration entry with proper multi-line indentation"""
    
    lines = [f"        '{code}': {{"]
    
    items = list(config.items())
    for idx, (key, value) in enumerate(items):
        is_last = (idx == len(items) - 1)
        comma = '' if is_last else ','
        
        if isinstance(value, dict):
            # Handle nested dictionary
            lines.append(f"        '{key}':")
            lines.append("          {")
            sub_items = list(value.items())
            for sub_idx, (sub_key, sub_value) in enumerate(sub_items):
                is_sub_last = (sub_idx == len(sub_items) - 1)
                sub_comma = '' if is_sub_last else ','
                lines.append(f"              '{sub_key}': {repr(sub_value)}{sub_comma}")
            lines.append("          },")
        else:
            lines.append(f"        '{key}': {repr(value)}{comma}")
    
    lines.append("        }")
    return '\n'.join(lines)

def _write_frame_entry_to_config(code, config):
    """Write a new frame entry into the Frame_Info dictionary in Columns_config.py"""
    with open('Columns_config.py', 'r') as f:
        content = f.read()
    
    # Format the entry with proper indentation
    formatted_entry = _format_frame_entry(code, config)
    
    # Format as a dictionary entry with comma, separator, and new entry
    entry_str = f",\n###########################################################################\n{formatted_entry}"
    
    # Find the last closing } of Frame_Info and insert before it
    closing_brace_index = content.rfind('}')
    if closing_brace_index != -1:
        # Remove trailing whitespace before the final }, so comma goes right after previous entry's }
        content = content[:closing_brace_index].rstrip() + entry_str + '\n' + content[closing_brace_index:]
    
    # Write back to file
    with open('Columns_config.py', 'w') as f:
        f.write(content)

def check_and_create_new_entries_in_column_config_file(column_section_name,story_height,no_of_stories,bending_axes,**kwargs):
    
    
    defaults={
        'bay_width': [],
        'column_no_of_ele': 2,
        'beam_no_of_ele': 4,
        'beam_section': {
            'common_and_exceptions': {
                'common': 'W27X84'
            }},
        'support': 'f',
        'load_comb_multipliers': [1, 0, 0, 1],
        'D_floor_intensity': 7.5 * 10 * KN,
        'D_roof_intensity': 133.5 *10* KN,
        'L_floor_intensity': 0 * kip / ft,
        'L_roof_intensity': 0 * kip / ft,
        'Wind_load_floor': 6.56 * 3 * KN,
        'Wind_load_roof': 4.45 * 3 * KN,
        'Wall_load': 0,
        'geometric_imperfection_ratio': 1 / 500,
        'Leaning_column': False,
        'Leaning_column_offset': 2,
        'Leaning_column_floor_load': 2,
        'Leaning_column_roof_load': 1
    }
    code = f"{column_section_name}_{bending_axes}_{no_of_stories}X{story_height}"
    if code not in Frame_Info:
        config = {}
        config['Frame_id']=code
        config['story_height']=[story_height]*no_of_stories
        config['column_section']={
                    'common_and_exceptions': {
                        'common': (column_section_name, bending_axes),
                    }}
        for key, value in defaults.items():
            config[key] = kwargs.get(key, value)
        Frame_Info[code]=config
        
        # Write to Columns_config.py as a proper dict entry within Frame_Info
        _write_frame_entry_to_config(code, config)
    return code     