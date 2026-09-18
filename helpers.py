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
from libdenavit.section.database.aisc import wide_flange_database as wide_flange_database
from libdenavit.interaction_diagram_2d import InteractionDiagram2d
import os
import json
from Columns_config import Frame_Info
from Analysis_config import Analysis_Info
from Materials_config import Material_Info
import ast
from pprint import pformat
import matplotlib as mpl
mpl.use("TkAgg")
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import math
import pandas as pd
import opsvis
from pathlib import Path
import imageio.v2 as imageio
import shutil
from Plots import line_plot
import numpy as np
import tempfile
from numbers import Real




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
            "wind_load_dirn_source": "unknown",
            "initial_out_of_straightness_dirn": None,
            "initial_out_of_straightness_dirn_source":"unknown"
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

def _save_frame_info_to_config():
    """Replace the existing Frame_Info dictionary in Columns_config.py."""

    file_path = "Columns_config.py"

    with open(file_path, "r", encoding="utf-8") as file:
        content = file.read()

    tree = ast.parse(content)

    frame_info_node = None

    # Find: Frame_Info = {...}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "Frame_Info":
                    frame_info_node = node.value
                    break

        if frame_info_node is not None:
            break

    if frame_info_node is None:
        raise ValueError(
            "Frame_Info dictionary was not found in Columns_config.py."
        )

    lines = content.splitlines(keepends=True)

    start_index = (
        sum(len(line) for line in lines[:frame_info_node.lineno - 1])
        + frame_info_node.col_offset
    )

    end_index = (
        sum(len(line) for line in lines[:frame_info_node.end_lineno - 1])
        + frame_info_node.end_col_offset
    )

    updated_frame_info = pformat(
        Frame_Info,
        sort_dicts=False,
        width=120
    )

    updated_content = (
        content[:start_index]
        + updated_frame_info
        + content[end_index:]
    )

    with open(file_path, "w", encoding="utf-8") as file:
        file.write(updated_content)

def  check_and_create_new_entries_in_column_config_file(
    column_section_name,
    story_height,
    no_of_stories,
    bending_axes,
    Material_type,
    Leaning_column,
    Floor_to_Roof_load_ratio,
    Leaning_Column_load_ratio,
    support,
    **kwargs
):
    defaults = {
        "bay_width": [],
        "column_no_of_ele": 4,
        "beam_no_of_ele": 4,
        "beam_section": {
            "common_and_exceptions": {
                "common": "W27X84"
            }
        },
        "support": support,
        "load_comb_multipliers": [1, 0, 0, 1],
        "floor_to_roof_load_ratio":Floor_to_Roof_load_ratio,
        "D_floor_intensity": 1000 * KN,
        "L_floor_intensity": 0 * kip / ft,
        "L_roof_intensity": 0 * kip / ft,
        "wind_load_same_for_all_h": False,
        "Base_Wind_load": 10 * KN,
        "Wall_load": 0,
        "geometric_imperfection_ratio": 1 / 500,
        "initial_out_of_straightness_ratio": 1/1000,
        "Leaning_column_offset": 2,
        "Leaning_Column_load_ratio":Leaning_Column_load_ratio,
        "floor_nodes_free": False

    }
    if Leaning_column:
        code = (
            f"LCLR_{Leaning_Column_load_ratio}_{column_section_name}_{bending_axes}_"
            f"{no_of_stories}_{story_height}_{Material_type}_support_{support}_FtR{Floor_to_Roof_load_ratio}"
        )
    else:
        code = (
            f"{column_section_name}_{bending_axes}_"
            f"{no_of_stories}_{story_height}_{Material_type}_support_{support}_FtR{Floor_to_Roof_load_ratio}"
        )

    config = {
        "Frame_id": code,
        "story_height": [story_height] * no_of_stories,
        "column_section": {
            "common_and_exceptions": {
                "common": (column_section_name, bending_axes)
            }
        },
        "Material_type": Material_type,
        "Leaning_column": Leaning_column,
    }

    for key, default_value in defaults.items():
        config[key] = kwargs.get(key, default_value)

    config["D_roof_intensity"]=config['D_floor_intensity']/Floor_to_Roof_load_ratio

    if Leaning_column:
        config['Leaning_column_floor_load']=config['D_floor_intensity']*Leaning_Column_load_ratio
        config['Leaning_column_roof_load']=config['Leaning_column_floor_load']/Floor_to_Roof_load_ratio
    else:
        config['Leaning_column_floor_load']=0
        config['Leaning_column_roof_load']=0


    # Case 1: The frame does not exist
    if code not in Frame_Info:
        Frame_Info[code] = config
        _save_frame_info_to_config()

        print(f"Created new configuration: {code}")

    # Case 2: The frame exists, but its values are different
    elif Frame_Info[code] != config:
        Frame_Info[code] = config
        _save_frame_info_to_config()

        print(f"Updated existing configuration: {code}")

    # Case 3: The frame exists and the values are unchanged
    else:
        print(f"Configuration already exists and is unchanged: {code}")

    return code   

def values_are_equal(
    existing_value,
    new_value,
    rel_tol=1e-9,
    abs_tol=1e-12
):
    """
    Recursively compare values at any nesting depth.

    Supported structures:
    - dictionaries
    - lists
    - tuples
    - numerical values
    - strings, booleans, and None

    Numerical values are compared using math.isclose().
    """

    # Compare dictionaries recursively
    if isinstance(existing_value, dict) and isinstance(new_value, dict):

        # Both dictionaries must contain exactly the same keys
        if existing_value.keys() != new_value.keys():
            return False

        return all(
            values_are_equal(
                existing_value[key],
                new_value[key],
                rel_tol=rel_tol,
                abs_tol=abs_tol
            )
            for key in existing_value
        )

    # Compare lists and tuples recursively
    if (
        isinstance(existing_value, (list, tuple))
        and isinstance(new_value, (list, tuple))
    ):

        # Both sequences must have the same number of elements
        if len(existing_value) != len(new_value):
            return False

        return all(
            values_are_equal(
                existing_item,
                new_item,
                rel_tol=rel_tol,
                abs_tol=abs_tol
            )
            for existing_item, new_item
            in zip(existing_value, new_value)
        )

    # Compare numerical values using tolerance
    # bool is excluded because bool is technically a subclass of int
    if (
        isinstance(existing_value, Real)
        and isinstance(new_value, Real)
        and not isinstance(existing_value, bool)
        and not isinstance(new_value, bool)
    ):
        return math.isclose(
            float(existing_value),
            float(new_value),
            rel_tol=rel_tol,
            abs_tol=abs_tol
        )

    # Compare strings, booleans, None, and other values normally
    return existing_value == new_value

def update_or_create_json(
    file_path,
    new_data,
    rel_tol=1e-5,
    abs_tol=1e-7,
    overwrite_on_conflict=True,
):
    """
    Create or update a JSON file safely.

    Rules
    -----
    1. Missing top-level keys are added.
    2. Equivalent existing values are left unchanged.
    3. Numerical values are compared with tolerance at any nesting depth.
    4. Conflicting top-level keys are:
       - rejected when overwrite_on_conflict=False;
       - replaced individually when overwrite_on_conflict=True.
    5. Unrelated existing keys are always preserved.
    """

    if not isinstance(new_data, dict) or not new_data:
        raise ValueError(
            "new_data must be a non-empty dictionary."
        )

    file_path = os.fspath(file_path)

    # Load existing JSON data
    if os.path.exists(file_path):
        try:
            with open(
                file_path,
                "r",
                encoding="utf-8",
            ) as json_file:
                existing_data = json.load(json_file)

        except json.JSONDecodeError as error:
            raise ValueError(
                f"'{file_path}' is empty or contains invalid JSON."
            ) from error

        if not isinstance(existing_data, dict):
            raise TypeError(
                f"The top-level data in '{file_path}' "
                "must be a dictionary."
            )

    else:
        existing_data = {}

    keys_to_add = []
    keys_to_overwrite = []
    equivalent_keys = []

    # First inspect everything without modifying existing_data
    for key, new_value in new_data.items():

        if key not in existing_data:
            keys_to_add.append(key)

        elif values_are_equal(
            existing_data[key],
            new_value,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
        ):
            equivalent_keys.append(key)

        else:
            keys_to_overwrite.append(key)

    # Reject conflicts when overwriting is disabled
    if keys_to_overwrite and not overwrite_on_conflict:
        conflict_text = ", ".join(
            repr(key) for key in keys_to_overwrite
        )

        raise ValueError(
            "Data discrepancy detected for key(s): "
            f"{conflict_text}.\n"
            "The JSON file was not changed.\n"
            "Set overwrite_on_conflict=True to replace only "
            "the conflicting key values."
        )

    # Nothing needs to change
    if not keys_to_add and not keys_to_overwrite:
        print(
            "All supplied keys already exist with equivalent data. "
            "No change required."
        )
        return

    # Preserve all existing keys
    data_to_write = existing_data.copy()

    # Add missing keys
    for key in keys_to_add:
        data_to_write[key] = new_data[key]

    # Replace only conflicting keys
    for key in keys_to_overwrite:
        data_to_write[key] = new_data[key]

    parent_folder = os.path.dirname(
        os.path.abspath(file_path)
    )
    os.makedirs(parent_folder, exist_ok=True)

    # Write safely using a temporary file
    temporary_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=parent_folder,
            delete=False,
            suffix=".json.tmp",
        ) as temporary_file:

            temporary_path = temporary_file.name

            json.dump(
                data_to_write,
                temporary_file,
                indent=2,
                separators=(",", ": "),
                ensure_ascii=False,
            )

        os.replace(
            temporary_path,
            file_path,
        )

    except Exception:
        if (
            temporary_path is not None
            and os.path.exists(temporary_path)
        ):
            os.remove(temporary_path)

        raise

    if equivalent_keys:
        print(
            "Equivalent keys left unchanged: "
            f"{equivalent_keys}"
        )

    if keys_to_add:
        print(
            f"Added {len(keys_to_add)} new key(s): "
            f"{keys_to_add}"
        )

    if keys_to_overwrite:
        print(
            f"Overwrote {len(keys_to_overwrite)} conflicting key(s): "
            f"{keys_to_overwrite}"
        )

    print(
        f"JSON file '{file_path}' updated successfully."
    )


def plot_interaction_diagrams_for_section_colormap(
    json_file_path,
    frame_name,
    no_of_stories_list=None,
    bending_axes=None,
    colormap="viridis",
    use_original=True,
    save=True,
    show=True
):
    """
    Plot interaction diagrams for all heights of a given section for a sanity check.

    """

    # save_folder=os.path.join("Column_Results",f"parametric{PARAMETER} vs max radial error")
    results_data = []
    Folder=Path(json_file_path)
    matched_results = []
    for file_path in Folder.iterdir():
        with open(file_path, "r", encoding="utf-8") as file:
            content = file.read()
        content=json.loads(content)
        analyses = [key for key in content.keys() if key != "data"]

        #### Read structural parameters from json file ###
        frame_id=content['data']['Frame_id'] 
        section_name=content['data']['column_section_name']
        bending_axis=content['data']["bending_axes"]
        storey_height=content['data']["storey_height"]
        no_of_stories=content['data']["no_of_stories"]


        matches_filter = ((section_name==frame_name) and
                        (no_of_stories_list is None or no_of_stories in no_of_stories_list)
                        and (bending_axes is None or bending_axis in bending_axes)
                        )

        if not matches_filter:
            continue
        matched_results.append(
            {
                "file_path": file_path,
                "content": content,
                "frame_id": frame_id,
                "section_name": frame_name,
                "bending_axis": bending_axis,
                "storey_height": storey_height,
                "no_of_stories": no_of_stories,
            })
        
        matched_results = sorted(
            matched_results,
            key=lambda item: item["storey_height"]
        )

    heights = np.array(
        [item["storey_height"] for item in matched_results],
        dtype=float
    )

    n_plots = len(analyses)
    n_cols = min(2, n_plots)
    n_rows = math.ceil(n_plots / n_cols)
    original_cmap = plt.get_cmap(colormap)

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        f"{colormap}_truncated",
        original_cmap(np.linspace(0, 0.5, 256))
    )

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(7 * n_cols, 5.5 * n_rows),
        sharex=True,
        sharey=True
    )

    axes = np.atleast_1d(axes).flatten()

    if np.isclose(heights.min(), heights.max()):
        norm = mpl.colors.Normalize(
            vmin=heights.min() - 0.5,
            vmax=heights.max() + 0.5
        )
    else:
        norm = mpl.colors.Normalize(
            vmin=heights.min(),
            vmax=heights.max()
        )


    for ax, analysis in zip(axes, analyses):

        for item in matched_results:
            content = item["content"]
            storey_height = item["storey_height"]
            section_name=item["section_name"]

            if analysis not in content:
                continue

            if use_original:
                alr_h_key = "Original_ALR_H"
                alr_v_key = "Original_ALR_V"
            else:
                alr_h_key = "ALR_H"
                alr_v_key = "ALR_V"

            if alr_h_key not in content[analysis] or alr_v_key not in content[analysis]:
                print(
                    f"Skipping {item['file_path'].name}: "
                    f"{analysis} does not contain {alr_h_key} and {alr_v_key}"
                )
                continue

            ALR_H = content[analysis][alr_h_key]
            ALR_V = content[analysis][alr_v_key]

            interaction_diagram = InteractionDiagram2d(ALR_H, ALR_V)

            # InteractionDiagram2d.plot uses the current axes
            plt.sca(ax)

            interaction_diagram.plot(
                color=cmap(norm(storey_height)),
                linewidth=1.3,
                alpha=0.9
            )

        ax.set_title(analysis)
        ax.set_xlabel("ALR_H")
        ax.set_ylabel("ALR_V")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)

    # Hide unused axes, if any
    for ax in axes[len(analyses):]:
        ax.set_visible(False)

    # ------------------------------------------------------------
    # 7. Add one shared colorbar for height
    # ------------------------------------------------------------
    scalar_mappable = mpl.cm.ScalarMappable(
        norm=norm,
        cmap=cmap
    )
    scalar_mappable.set_array([])

    # Leave space on the right side for the colorbar
    fig.tight_layout(rect=[0, 0, 0.86, 0.95])

    # Create a separate axis for the colorbar
    cbar_ax = fig.add_axes([0.95, 0.15, 0.015, 0.50])

    colorbar = fig.colorbar(
        scalar_mappable,
        cax=cbar_ax
    )

    colorbar.set_label("Storey height, h (m)")

    if len(np.unique(heights)) > 6:
        colorbar.set_ticks(
            np.linspace(heights.min(), heights.max(), 6)
        )

    # ------------------------------------------------------------
    # 8. Overall title
    # ------------------------------------------------------------
    title = f"{section_name}"

    if bending_axes is not None and len(bending_axes) == 1:
        title += f"_{bending_axes[0]}"

    if no_of_stories_list is not None and len(no_of_stories_list) == 1:
        title += f"_{no_of_stories_list[0]} story"

    fig.suptitle(title, fontsize=14)

    if save:
        save_folder = os.path.join(
            "Column_Results",
            "Check_Interaction_Diagrams",
            f"{frame_name}"
        )
        os.makedirs(save_folder, exist_ok=True)

        filename = f"{section_name}_interaction_diagrams_all_heights"

        if bending_axes is not None:
            filename += "_" + "_".join(str(axis) for axis in bending_axes)

        if no_of_stories_list is not None:
            filename += "_" + "_".join(
                f"{story}story" for story in no_of_stories_list
            )

        filename += ".png"

        save_path = os.path.join(save_folder, filename)

        fig.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight"
        )

        print(f"Saved to: {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig, axes

def plot_interaction_diagrams_for_all_analyses(
    json_file_path,
    use_original=True,
    save=True,
    show=True,
):
    """
    Read every JSON file inside a folder and plot all analysis interaction
    diagrams from each JSON file on the same figure.

    One figure is created for each JSON file.

    Parameters
    ----------
    json_file_path : str or Path
        Path to the folder containing the JSON files.

    use_original : bool, default=True
        If True, plot Original_ALR_H and Original_ALR_V.
        If False, plot ALR_H and ALR_V.

    save : bool, default=True
        Save each generated figure.

    show : bool, default=True
        Display each generated figure.
    """

    Folder = Path(json_file_path)

    if not Folder.exists():
        raise FileNotFoundError(f"The specified path does not exist: {Folder}")

    if not Folder.is_dir():
        raise NotADirectoryError(f"The specified path is not a directory: {Folder}")

    json_files = sorted(
        file_path
        for file_path in Folder.iterdir()
        if file_path.is_file() and file_path.suffix.lower() == ".json"
    )

    if not json_files:
        print(f"No JSON files were found in: {Folder}")
        return []

    if use_original:
        alr_h_key = "Original_ALR_H"
        alr_v_key = "Original_ALR_V"
    else:
        alr_h_key = "ALR_H"
        alr_v_key = "ALR_V"

    results = []

    for file_number, file_path in enumerate(json_files, start=1):
        try:
            with open(file_path, "r", encoding="utf-8") as file:
                content = json.load(file)

        except json.JSONDecodeError as error:
            print(f"Skipping invalid JSON file {file_path.name}: {error}")
            continue

        except OSError as error:
            print(f"Could not read {file_path.name}: {error}")
            continue

        data = content.get("data", {})

        frame_id = data.get("frame_id", file_path.stem)
        section_name = data.get("column_section_name", "Unknown section")
        bending_axis = data.get("bending_axes", data.get("bending_axis", "Unknown axis"))
        storey_height = data.get("storey_height", "Unknown height")
        no_of_stories = data.get("no_of_stories", "Unknown number of stories")
        material = data.get("Material", data.get("material", None))

        analyses = [
            key
            for key, value in content.items()
            if key != "data" and isinstance(value, dict)
        ]

        if not analyses:
            print(f"Skipping {file_path.name}: no analyses were found.")
            continue

        fig, ax = plt.subplots(figsize=(8, 6))
        plotted_analyses = []
        line_styles = ["-", "--", "-.", ":"]

        for analysis_index, analysis in enumerate(analyses):
            analysis_data = content[analysis]

            if alr_h_key not in analysis_data or alr_v_key not in analysis_data:
                print(
                    f"Skipping analysis '{analysis}' in {file_path.name}: "
                    f"missing {alr_h_key} or {alr_v_key}."
                )
                continue

            ALR_H = analysis_data[alr_h_key]
            ALR_V = analysis_data[alr_v_key]

            if ALR_H is None or ALR_V is None:
                print(
                    f"Skipping analysis '{analysis}' in {file_path.name}: "
                    "interaction data are None."
                )
                continue

            if len(ALR_H) == 0 or len(ALR_V) == 0:
                print(
                    f"Skipping analysis '{analysis}' in {file_path.name}: "
                    "interaction data are empty."
                )
                continue

            if len(ALR_H) != len(ALR_V):
                print(
                    f"Skipping analysis '{analysis}' in {file_path.name}: "
                    "ALR_H and ALR_V have different lengths."
                )
                continue

            interaction_diagram = InteractionDiagram2d(ALR_H, ALR_V)

            plt.sca(ax)
            interaction_diagram.plot(
                label=analysis,
                linestyle=line_styles[analysis_index % len(line_styles)],
                linewidth=1.8,
                alpha=0.9,
            )

            plotted_analyses.append(analysis)

        if not plotted_analyses:
            print(f"Skipping {file_path.name}: no valid interaction curves were found.")
            plt.close(fig)
            continue

        ax.set_xlabel("ALR_H")
        ax.set_ylabel("ALR_V")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        ax.legend(title="Analysis", loc="best")

        if no_of_stories == 1:
            story_label = "1 story"
        elif isinstance(no_of_stories, (int, float)):
            story_label = f"{no_of_stories} stories"
        else:
            story_label = str(no_of_stories)

        title_lines = [
            f"{section_name}, axis = {bending_axis}",
            f"{story_label}, h = {storey_height} m",
        ]

        if material is not None:
            title_lines.append(f"Material = {material}")

        ax.set_title("\n".join(title_lines))
        fig.tight_layout()

        save_path = None

        if save:
            save_folder = os.path.join("Column_Results/Analysis_History", f"{frame_id}")
            os.makedirs(save_folder, exist_ok=True)

            filename = f"{frame_id}_interaction_diagrams_all_analyses.png"
            save_path = os.path.join(save_folder, filename)

            fig.savefig(save_path, dpi=300, bbox_inches="tight")

            print(f"[{file_number}/{len(json_files)}] Saved to: {save_path}")

        if show:
            plt.show()
        else:
            plt.close(fig)

def plot_theta_vs_del2_over_del1_for_all_analyses(json_file_path,cap=3, save=True, show=True):
    """
    Read every JSON file inside a folder and plot theta vs del2_over_del1
    for all analyses from each JSON file on the same figure.

    One figure is created for each JSON file.
    """

    Folder = Path(json_file_path)

    if not Folder.exists():
        raise FileNotFoundError(f"The specified path does not exist: {Folder}")

    if not Folder.is_dir():
        raise NotADirectoryError(f"The specified path is not a directory: {Folder}")

    json_files = sorted(file_path for file_path in Folder.iterdir() if file_path.is_file() and file_path.suffix.lower() == ".json")

    if not json_files:
        print(f"No JSON files were found in: {Folder}")
        return []

    results = []

    for file_number, file_path in enumerate(json_files, start=1):

        try:
            with open(file_path, "r", encoding="utf-8") as file:
                content = json.load(file)

        except json.JSONDecodeError as error:
            print(f"Skipping invalid JSON file {file_path.name}: {error}")
            continue

        except OSError as error:
            print(f"Could not read {file_path.name}: {error}")
            continue

        data = content.get("data", {})

        frame_id = data.get("Frame_id", file_path.stem)
        section_name = data.get("column_section_name", "Unknown section")
        bending_axis = data.get("bending_axes", data.get("bending_axis", "Unknown axis"))
        storey_height = data.get("storey_height", "Unknown height")
        no_of_stories = data.get("no_of_stories", "Unknown number of stories")
        material = data.get("Material", data.get("material", None))

        analyses = [key for key, value in content.items() if key != "data" and isinstance(value, dict)]

        if not analyses:
            print(f"Skipping {file_path.name}: no analyses were found.")
            continue

        fig, ax = plt.subplots(figsize=(8, 6))
        plotted_analyses = []
        line_styles = ["-", "--", "-.", ":"]

        for analysis_index, analysis in enumerate(analyses):

            analysis_data = content[analysis]

            if "Theta" not in analysis_data or "del2_over_del1" not in analysis_data:
                print(f"Skipping analysis '{analysis}' in {file_path.name}: missing theta or del2_over_del1.")
                continue

            theta = analysis_data['del2_over_del1_results']["Theta"]
            del2_over_del1 = analysis_data['del2_over_del1_results']["del2_over_del1"]

            del2_over_del1_capped=[min(x,cap) for x in del2_over_del1]

            if theta is None or del2_over_del1 is None:
                print(f"Skipping analysis '{analysis}' in {file_path.name}: theta or del2_over_del1 is None.")
                continue

            if len(theta) == 0 or len(del2_over_del1) == 0:
                print(f"Skipping analysis '{analysis}' in {file_path.name}: theta or del2_over_del1 is empty.")
                continue

            if len(theta) != len(del2_over_del1):
                print(f"Skipping analysis '{analysis}' in {file_path.name}: theta and del2_over_del1 have different lengths.")
                continue

            ax.plot(theta, del2_over_del1_capped, label=analysis, linestyle=line_styles[analysis_index % len(line_styles)], linewidth=1.8, alpha=0.9,marker='o', markersize=2.5)

            plotted_analyses.append(analysis)

        if not plotted_analyses:
            print(f"Skipping {file_path.name}: no valid theta vs del2_over_del1 curves were found.")
            plt.close(fig)
            continue

        ax.set_xlabel(r"$\theta$")
        ax.set_ylabel(r"$\Delta_2/\Delta_1$")
        ax.grid(True, alpha=0.3)
        ax.legend(title="Analysis", loc="best")

        if no_of_stories == 1:
            story_label = "1 story"
        elif isinstance(no_of_stories, (int, float)):
            story_label = f"{no_of_stories} stories"
        else:
            story_label = str(no_of_stories)

        title_lines = [f"{section_name}, axis = {bending_axis}", f"{story_label}, h = {storey_height} m"]

        if material is not None:
            title_lines.append(f"Material = {material}")

        ax.set_title("\n".join(title_lines))
        fig.tight_layout()

        save_path = None

        if save:
            save_folder = os.path.join("Column_Results", "Analysis_History", f"{frame_id}")
            os.makedirs(save_folder, exist_ok=True)

            filename = f"{frame_id}_theta_vs_del2_over_del1_all_analyses.png"
            save_path = os.path.join(save_folder, filename)

            fig.savefig(save_path, dpi=300, bbox_inches="tight")
            print(f"[{file_number}/{len(json_files)}] Saved to: {save_path}")

        results.append({"frame_id": frame_id, "file_path": str(file_path), "analyses": plotted_analyses, "save_path": save_path})

        if show:
            plt.show()
        else:
            plt.close(fig)

    return results




def my_scatter_plot(data_df, x_col, y_col, hover_data=None, xlabel=None, ylabel=None, title='Interactive Scatter Plot'):

    if xlabel is None:
        xlabel = x_col
    if ylabel is None:
        ylabel = y_col
    if hover_data is None:
        hover_data = True
    
    fig = px.scatter(data_df, 
                     x=x_col, 
                     y=y_col,
                     hover_data=hover_data,
                     title=title,
                     labels={x_col: xlabel, y_col: ylabel})
    fig.update_traces(
        marker=dict(
            size=3,
            opacity=0.9,
            line=dict(color="black", width=0.7)
        )
    )

    fig.update_xaxes(
        showgrid=True,
        gridcolor="rgba(0,0,0,0.12)",
        gridwidth=1,
        showline=True,
        linewidth=1.3,
        linecolor="black",
        mirror=True,
        ticks="outside",
        ticklen=5,
        tickwidth=1,
        tickcolor="black",
        zeroline=True,
        zerolinecolor="black",
        zerolinewidth=2.2
    )

    fig.update_yaxes(
        showgrid=True,
        gridcolor="rgba(0,0,0,0.12)",
        gridwidth=1,
        showline=True,
        linewidth=1.3,
        linecolor="black",
        mirror=True,
        ticks="outside",
        ticklen=5,
        tickwidth=1,
        tickcolor="black",
        zeroline=True,
        zerolinecolor="black",
        zerolinewidth=2.2
    )

    fig.update_layout(
        width=800,
        height=550,
        plot_bgcolor="white",
        paper_bgcolor="white",
        font=dict(family="Arial", size=15, color="black"),
        title=dict(x=0.5, xanchor="center"),
        margin=dict(l=75, r=20, t=60, b=65),
    )



    return fig

def structural_parameter_vs_radial_error(json_file_path, analysis1, analysis2,PARAMETER='slenderness_ratio', plot_interaction_diagram=False,
                                      frame_names_list=None,
                                      storey_heights_list=None,
                                      no_of_stories_list=None,
                                      bending_axes=None, **kwargs):

    save_folder=os.path.join("Column_Results",f"parametric{PARAMETER} vs max radial error")
    results_data = []
 
    Folder=Path(json_file_path)
    for file_path in Folder.iterdir():
        with open(file_path, "r", encoding="utf-8") as file:
            content = file.read()
        content=json.loads(content)

        #### Read structural parameters from json file ###
        frame_id=content['data']['Frame_id'] 
        frame_name=content['data']['column_section_name']
        bending_axis=content['data']["bending_axes"]
        storey_height=content['data']["storey_height"]
        no_of_stories=content['data']["no_of_stories"]
        support=content['data']["support"]


        #### Filter for the lists if provided ####
        matches_filter = (
                        (frame_names_list is None or frame_name in frame_names_list)
                        and (storey_heights_list is None or storey_height in storey_heights_list)
                        and (no_of_stories_list is None or no_of_stories in no_of_stories_list)
                        and (bending_axes is None or bending_axis in bending_axes)
                        )

        if not matches_filter:
            continue

        ##### Get section properties from aisc database ####
        r=(wide_flange_database[frame_name]['rx'] if bending_axis=='x' else wide_flange_database[frame_name]['ry'])* inch


        # Extract ALR values for both analyses
        analysis1_ALR_H = content[analysis1]['Original_ALR_H']
        analysis1_ALR_V = content[analysis1]['Original_ALR_V']
        analysis2_ALR_H = content[analysis2]['Original_ALR_H']
        analysis2_ALR_V = content[analysis2]['Original_ALR_V']
    
        analysis1_interaction_diagram = InteractionDiagram2d(analysis1_ALR_H, analysis1_ALR_V)
        analysis2_interaction_diagram = InteractionDiagram2d(analysis2_ALR_H, analysis2_ALR_V)
        


        theta_list = np.linspace(0, np.pi/2, 91) 
        errors = analysis1_interaction_diagram.compare_two(analysis2_interaction_diagram, theta_list, degrees=False)

        max_positive_error=np.max(errors)
        id_max_positive=np.argmax(errors)
        min_negative_error=np.min(errors)
        id_min_negative=np.argmin(errors)
        theta_min_error = theta_list[id_min_negative]
        theta_max_error = theta_list[id_max_positive]

        d_min=analysis2_interaction_diagram.radial_distance( theta_min_error, degrees=False)
        d_max=analysis1_interaction_diagram.radial_distance(theta_max_error,degrees=False)


        if plot_interaction_diagram:

            fig_interaction, ax = plt.subplots(figsize=(7, 6))
            plt.sca(ax)

            analysis1_interaction_diagram.plot('-r', label=analysis1)
            analysis2_interaction_diagram.plot('-b', label=analysis2)

            x_min=[0,d_min*np.cos(theta_min_error)]
            y_min=[0,d_min*np.sin(theta_min_error)]

            x_max=[0,d_max*np.cos(theta_max_error)]
            y_max=[0,d_max*np.sin(theta_max_error)]

            ax.plot(x_min, y_min, color='green', linestyle='--', linewidth=2, label=f'Max Unconservative Error ({min_negative_error:.4f})')

            ax.plot(x_max, y_max, color='magenta', linestyle='-.', linewidth=2, label=f'Max Conservative Error ({max_positive_error:.4f})')

            plt.legend()
            plt.title(f'Interaction Diagram for Frame {frame_id}')
            plt.xlabel('ALR_H')
            plt.ylabel('ALR_V')
            plt.grid()
            interaction_diagram_save_path=os.path.join("Column_Results/Analysis_History", f"{frame_id}")
            os.makedirs(interaction_diagram_save_path, exist_ok=True)
            plt.savefig(os.path.join(interaction_diagram_save_path, f'interaction_diagram_{frame_id}_{analysis1}_vs_{analysis2}.png'))
            plt.close()

        if PARAMETER.lower()=='slenderness_ratio':
            parameter = storey_height / r
        
        # Store data as dictionary for DataFrame
        results_data.append({
            'frame': frame_id,
            'section': frame_name,
            'axis': bending_axis,
            'stories': no_of_stories,
            'height': storey_height,
            'support': support,
            'radius_of_gyration': r,
            f'{PARAMETER}': parameter,
            'max_positive_error': max_positive_error,
            'min_negative_error':min_negative_error
        })
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results_data)
    
    # Define hover data with formatting
    hover_data_config = {
        'section': True, 
        'axis': True, 
        'stories': True, 
        'height': ':.3f',
        'support': support,
        'radius_of_gyration': ':.3f',
        f'{PARAMETER}': ':.3f', 
        'max_positive_error': ':.6f',
        'min_negative_error':':.6f'

    }
    
    # Create interactive Plotly scatter plot
    fig = my_scatter_plot(results_df, 
                       x_col=f'{PARAMETER}',
                       y_col='min_negative_error',
                       hover_data=hover_data_config,
                       xlabel=f'{PARAMETER}',
                       ylabel='Minimum Radial Error',
                       title=f'{PARAMETER} vs Radial Error ({analysis1} > {analysis2})')
    
    os.makedirs(save_folder, exist_ok=True)
    
    # Define common base filename
    base_filename = f"{PARAMETER}_vs_min_radial_error_{analysis1}_{analysis2}"
    
    # 1. Save as interactive HTML 
    html_path = os.path.join(save_folder, f"{base_filename}.html")
    fig.write_html(html_path)
    print(f"Interactive plot saved to: {html_path}")
    
    png_path = os.path.join(save_folder, f"{base_filename}.png")
    
    fig.write_image(png_path, width=1000, height=600) 
    print(f"Static PNG image saved to: {png_path}")
    fig.show()
    
    return results_df, fig

def structural_parameter_vs_radial_error_multiple(json_file_path, analysis_types, PARAMETER='slenderness_ratio', plot_interaction_diagram=False,
                                      frame_names_list=None,
                                      storey_heights_list=None,
                                      no_of_stories_list=None,
                                      bending_axes_list=None,
                                      support_list=None, **kwargs):

    if len(analysis_types) < 2:
        print(analysis_types)
        raise ValueError("At least two analysis types must be provided.")

    analysis1 = analysis_types[0]

    save_folder=os.path.join("Column_Results",f"parametric{PARAMETER} vs max radial error")
    results_data = []
 
    Folder=Path(json_file_path)
    for file_path in Folder.iterdir():
        with open(file_path, "r", encoding="utf-8") as file:
            content = file.read()
        content=json.loads(content)

        #### Read structural parameters from json file ###
        frame_id=content['data']['Frame_id'] 
        frame_name=content['data']['column_section_name']
        bending_axis=content['data']["bending_axes"]
        storey_height=content['data']["storey_height"][0]
        no_of_stories=content['data']["no_of_stories"]
        support=content['data']["support"]


        #### Filter for the lists if provided ####
        matches_filter = (
                        (frame_names_list is None or frame_name in frame_names_list)
                        and (storey_heights_list is None or storey_height in storey_heights_list)
                        and (no_of_stories_list is None or no_of_stories in no_of_stories_list)
                        and (bending_axes_list is None or bending_axis in bending_axes_list)
                        and (support_list is None or support in support_list)
                        )

        if not matches_filter:
            continue

        ##### Get section properties from aisc database ####
        r=(wide_flange_database[frame_name]['rx'] if bending_axis=='x' else wide_flange_database[frame_name]['ry'])* inch


        # Extract ALR values for reference analysis
        analysis1_ALR_H = content[analysis1]['Original_ALR_H']
        analysis1_ALR_V = content[analysis1]['Original_ALR_V']
    
        analysis1_interaction_diagram = InteractionDiagram2d(analysis1_ALR_H, analysis1_ALR_V)


        #####################################################################
        # Compare analysis1 with every other analysis provided
        #####################################################################

        for analysis2 in analysis_types[1:]:

            analysis2_ALR_H = content[analysis2]['Original_ALR_H']
            analysis2_ALR_V = content[analysis2]['Original_ALR_V']

            analysis2_interaction_diagram = InteractionDiagram2d(analysis2_ALR_H, analysis2_ALR_V)
            
            theta_list = np.linspace(0, 90, 1000) 
            errors = analysis1_interaction_diagram.compare_two(analysis2_interaction_diagram, theta_list, degrees=True)

            # plt.figure(figsize=(7, 6))
            # plt.plot(theta_list, errors)
            # plt.xlabel('Theta (radians)')
            # plt.ylabel('Radial Error')
            # plt.title(f'Radial Error vs Theta ({analysis1} vs {analysis2})')
            # plt.grid()
            # plt.show()

            # input()
            

            max_positive_error=np.max(errors)
            id_max_positive=np.argmax(errors)
            min_negative_error=np.min(errors)
            id_min_negative=np.argmin(errors)
            theta_min_error = theta_list[id_min_negative]
            theta_max_error = theta_list[id_max_positive]

            d_min=analysis2_interaction_diagram.radial_distance(theta_min_error, degrees=True)
            d_max=analysis1_interaction_diagram.radial_distance(theta_max_error,degrees=True)


            if plot_interaction_diagram:

                fig_interaction, ax = plt.subplots(figsize=(7, 6))
                plt.sca(ax)

                analysis1_interaction_diagram.plot('-r', label=analysis1)
                analysis2_interaction_diagram.plot('-b', label=analysis2)

                x_min=[0,d_min*np.cos(np.radians(theta_min_error))]
                y_min=[0,d_min*np.sin(np.radians(theta_min_error))]

                x_max=[0,d_max*np.cos(np.radians(theta_max_error))]
                y_max=[0,d_max*np.sin(np.radians(theta_max_error))]

                ax.plot(x_min, y_min, color='green', linestyle='--', linewidth=2, label=f'Max Unconservative Error ({min_negative_error:.4f})')

                ax.plot(x_max, y_max, color='magenta', linestyle='-.', linewidth=2, label=f'Max Conservative Error ({max_positive_error:.4f})')

                plt.legend()
                plt.title(f'Interaction Diagram for Frame {frame_id}')
                plt.xlabel('ALR_H')
                plt.ylabel('ALR_V')
                plt.grid()
                interaction_diagram_save_path=os.path.join("Column_Results/Analysis_History", f"{frame_id}")
                os.makedirs(interaction_diagram_save_path, exist_ok=True)
                plt.savefig(os.path.join(interaction_diagram_save_path, f'interaction_diagram_{frame_id}_{analysis1}_vs_{analysis2}.png'))
                plt.close()

            if PARAMETER.lower()=='slenderness_ratio':
                parameter = storey_height / r
            
            # Store data as dictionary for DataFrame
            results_data.append({
                'frame': frame_id,
                'section': frame_name,
                'axis': bending_axis,
                'stories': no_of_stories,
                'height': storey_height,
                'support': support,
                'radius_of_gyration': r,
                f'{PARAMETER}': parameter,
                'comparison': f'{analysis1} vs {analysis2}',
                'max_positive_error': max_positive_error,
                'min_negative_error':min_negative_error
            })
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results_data)
    
    # Define hover data with formatting
    hover_data_config = {
        'section': True, 
        'axis': True, 
        'stories': True, 
        'height': ':.3f',
        'support': True,
        'radius_of_gyration': ':.3f',
        f'{PARAMETER}': ':.3f',
        'comparison': True,
        'max_positive_error': ':.6f',
        'min_negative_error':':.6f'
    }


    #####################################################################
    # Create one scatter plot for each comparison and combine them
    #####################################################################

    fig = None

    colors = px.colors.qualitative.Plotly

    for i, analysis2 in enumerate(analysis_types[1:]):

        comparison_name = f'{analysis1} vs {analysis2}'

        comparison_df = results_df[
            results_df['comparison'] == comparison_name
        ]

        temp_fig = my_scatter_plot(
            comparison_df, 
            x_col=f'{PARAMETER}',
            y_col='min_negative_error',
            hover_data=hover_data_config,
            xlabel=f'{PARAMETER}',
            ylabel='Minimum Radial Error',
            title=f'{PARAMETER} vs Radial Error ({analysis1} > Other Analyses)'
        )

        # Give this comparison its own color and legend name
        for trace in temp_fig.data:
            trace.name = comparison_name
            trace.showlegend = True
            trace.marker.color = colors[i % len(colors)]

        # First comparison becomes the main figure
        if fig is None:
            fig = temp_fig

        # Add subsequent comparisons to the same figure
        else:
            for trace in temp_fig.data:
                fig.add_trace(trace)


    os.makedirs(save_folder, exist_ok=True)
    
    # Define common base filename
    base_filename = f"{PARAMETER}_vs_min_radial_error_{analysis1}_vs_{'_'.join(analysis_types[1:])}"
    
    # 1. Save as interactive HTML 
    html_path = os.path.join(save_folder, f"{base_filename}.html")
    fig.write_html(html_path)
    print(f"Interactive plot saved to: {html_path}")
    
    png_path = os.path.join(save_folder, f"{base_filename}.png")
    
    fig.write_image(png_path, width=1000, height=600) 
    print(f"Static PNG image saved to: {png_path}")
    fig.show()
    
    return results_df, fig


def compute_del2_over_del1(json_file_path,theta_list):

    from MFColumn2D import MFColumn_2D
    """
    Read every JSON file inside a folder and compute del2_over_del1 from the original ALR_H and ALR_V values.


    Parameters
    ----------
    json_file_path : str or Path
        Path to the folder containing the JSON files.

    """

    Folder = Path(json_file_path)

    if not Folder.exists():
        raise FileNotFoundError(f"The specified path does not exist: {Folder}")

    if not Folder.is_dir():
        raise NotADirectoryError(f"The specified path is not a directory: {Folder}")

    json_files = sorted(
        file_path
        for file_path in Folder.iterdir()
        if file_path.is_file() and file_path.suffix.lower() == ".json"
    )

    if not json_files:
        print(f"No JSON files were found in: {Folder}")
        return []

    alr_h_key = "Original_ALR_H"
    alr_v_key = "Original_ALR_V"




    for file_number, file_path in enumerate(json_files, start=1):
        try:
            with open(file_path, "r", encoding="utf-8") as file:
                content = json.load(file)

        except json.JSONDecodeError as error:
            print(f"Skipping invalid JSON file {file_path.name}: {error}")
            continue

        except OSError as error:
            print(f"Could not read {file_path.name}: {error}")
            continue

        data = content.get("data", {})

        frame_id = data.get("Frame_id", file_path.stem)
        width_of_bay = data.get("width_of_bay", [])
        storey_height = data.get("storey_height", [])
        no_of_elements_column = data.get("no_of_elements_column", None)
        no_of_elements_beam = data.get("no_of_elements_beam", None)

        beam_section = data.get("beam_section", None)
        column_section = data.get("column_section", None)
        column_section_name = data.get("column_section_name", "Unknown section")
        bending_axis = data.get("bending_axes", data.get("bending_axis", "Unknown axis"))
        no_of_stories = data.get("no_of_stories", None)

        support = data.get("support", None)

        D_floor_intensity = data.get("D_floor_intensity", None)
        D_roof_intensity = data.get("D_roof_intensity", None)
        L_floor_intensity = data.get("L_floor_intensity", None)
        L_roof_intensity = data.get("L_roof_intensity", None)
        Base_Wind_load = data.get("Base_Wind_load", None)
        Wall_load = data.get("Wall_load", None)

        load_combination_multipliers = data.get("load_combination_multipliers", None)

        Steel_Grade = data.get("Steel_Grade", None)

        Residual_Stress = data.get("Residual_Stress", None)
        Elastic_analysis = data.get("Elastic_analysis", None)
        Second_order_effects = data.get("Second_order_effects", None)
        stiffness_reduction = data.get("stiffness_reduction", None)
        strength_reduction = data.get("strength_reduction", None)
        Notional_load = data.get("Notional_load", None)
        Geometric_Imperfection = data.get("Geometric_Imperfection", None)
        geometric_imperfection_ratio = data.get("geometric_imperfection_ratio", None)
        initial_out_of_straightness_ratio=data.get("initial_out_of_straightness_ratio", None)
        initial_out_of_straightness_dirn=data.get("initial_out_of_straightness_dirn", None)
        nip = data.get("nip", None)
        mat_type = data.get("mat_type", None)
        wind_load_dirn = data.get("wind_load_dirn", None)

        Leaning_column = data.get("Leaning_column", None)
        Leaning_column_offset = data.get("Leaning_column_offset", None)
        Leaning_column_floor_load = data.get("Leaning_column_floor_load", None)
        Leaning_column_roof_load = data.get("Leaning_column_roof_load", None)

        floor_nodes_free = data.get("floor_nodes_free", None)
        wind_load_same_for_all_h = data.get("wind_load_same_for_all_h", None)

        analyses = [
            key
            for key, value in content.items()
            if key != "data" and isinstance(value, dict)
        ]

        if not analyses:
            print(f"Skipping {file_path.name}: no analyses were found.")
            continue

        for analysis_index, analysis in enumerate(analyses):
            analysis_data = content[analysis]


            if alr_h_key not in analysis_data or alr_v_key not in analysis_data:
                print(
                    f"Skipping analysis '{analysis}' in {file_path.name}: "
                    f"missing {alr_h_key} or {alr_v_key}."
                )
                continue

            ALR_H = analysis_data[alr_h_key]
            ALR_V = analysis_data[alr_v_key]

            if ALR_H is None or ALR_V is None:
                print(
                    f"Skipping analysis '{analysis}' in {file_path.name}: "
                    "interaction data are None."
                )
                continue

            if len(ALR_H) == 0 or len(ALR_V) == 0:
                print(
                    f"Skipping analysis '{analysis}' in {file_path.name}: "
                    "interaction data are empty."
                )
                continue

            if len(ALR_H) != len(ALR_V):
                print(
                    f"Skipping analysis '{analysis}' in {file_path.name}: "
                    "ALR_H and ALR_V have different lengths."
                )
                continue

            interaction_curve = InteractionDiagram2d(ALR_H, ALR_V)

            results={
                    'Theta':[],
                    'ALR_H':[],
                    'ALR_V':[],
                    'del2_over_del1':[]
                    }
            
            for theta in theta_list:
                pathX, pathY = ((0, 0), (0, 1)) if theta == 90 else ((0, 1), (0, math.tan(math.radians(theta))))
                pt = interaction_curve.find_intersection(pathX, pathY)
                vertical_load_scale = pt[1] if pt is not None else None
                if vertical_load_scale==0:
                    ult_lat_load_for_V0=pt[0]
                    no_of_digits = len(str(int(abs(ult_lat_load_for_V0)))) if ult_lat_load_for_V0!=0 else 1
                lateral_load_scale = ult_lat_load_for_V0*(10**(-no_of_digits - 2)) if pt[0]<10e-6 else pt[0] if pt is not None else None
                print(Steel_Grade)
                Material_dict=Material_Info[Steel_Grade]
                Material_details=convert_dict_items_to_class_attributes(Material_dict)
                Steel=Steel_Material(mat_tag=1,E=Material_details.E,Fy=Material_details.Fy)

                Frame = MFColumn_2D(width_of_bay, storey_height, no_of_elements_column, no_of_elements_beam,
                                    beam_section=beam_section,
                                    column_section=column_section,
                                    support=support,
                                    D_floor_intensity=D_floor_intensity,
                                    D_roof_intensity=D_roof_intensity,
                                    L_floor_intensity=L_floor_intensity,
                                    L_roof_intensity=L_roof_intensity,
                                    Base_Wind_load=Base_Wind_load,
                                    Wall_load=Wall_load,
                                    load_combination_multipliers=load_combination_multipliers,
                                    Frame_id=frame_id,
                                    Material_obj=Steel,
                                    Steel_Grade=Steel_Grade,
                                    Residual_Stress=Residual_Stress,
                                    Elastic_analysis=Elastic_analysis,
                                    Second_order_effects=Second_order_effects,
                                    stiffness_reduction=stiffness_reduction,
                                    strength_reduction=strength_reduction,
                                    Notional_load=Notional_load,
                                    Geometric_Imperfection=Geometric_Imperfection,
                                    geometric_imperfection_ratio=geometric_imperfection_ratio,
                                    initial_out_of_straightness_ratio=initial_out_of_straightness_ratio,
                                    initial_out_of_straightness_dirn=initial_out_of_straightness_dirn,
                                    nip=nip,
                                    mat_type=mat_type,
                                    wind_load_dirn=wind_load_dirn,
                                    Leaning_column=Leaning_column,
                                    Leaning_column_offset=Leaning_column_offset,
                                    Leaning_column_floor_load=Leaning_column_floor_load,
                                    Leaning_column_roof_load=Leaning_column_roof_load,
                                    floor_nodes_free=floor_nodes_free,
                                    wind_load_same_for_all_h=wind_load_same_for_all_h)



                del2_over_del1=Frame.get_del2_over_del1(vertical_load_scale=vertical_load_scale, lateral_load_scale=lateral_load_scale)
                results['Theta'].append(theta)
                results['ALR_H'].append(pt[0] if pt is not None else None)
                results['ALR_V'].append(pt[1] if pt is not None else None)
                results['del2_over_del1'].append(del2_over_del1)

            updated_analysis_data = analysis_data.copy()
            updated_analysis_data["del2_over_del1_results"] = results

            Results = {analysis: updated_analysis_data}
            update_or_create_json(file_path, Results)

            
def plot_error_vs_del2_over_del1(json_file_path, analysis_types,
                                 cap=3,
                                frame_names_list=None,
                                storey_heights_list=None,
                                no_of_stories_list=None,
                                bending_axes_list=None,
                                support_list=None, **kwargs):
    
    if len(analysis_types) < 2:
        print(analysis_types)
        raise ValueError("At least two analysis types must be provided.")

    analysis1 = analysis_types[0]

    save_folder = os.path.join("Column_Results", "error_vs_del2_over_del1")
    results_data = []

    Folder = Path(json_file_path)

    for file_path in Folder.iterdir():

        if not file_path.is_file() or file_path.suffix.lower() != ".json":
            continue

        with open(file_path, "r", encoding="utf-8") as file:
            content = json.load(file)

        #### Read structural parameters from json file ###
        frame_id = content['data']['Frame_id']
        frame_name = content['data']['column_section_name']
        bending_axis = content['data']["bending_axes"]
        storey_height = content['data']["storey_height"][0]
        no_of_stories = content['data']["no_of_stories"]
        support = content['data']["support"]

        #### Filter for the lists if provided ####
        matches_filter = (
            (frame_names_list is None or frame_name in frame_names_list)
            and (storey_heights_list is None or storey_height in storey_heights_list)
            and (no_of_stories_list is None or no_of_stories in no_of_stories_list)
            and (bending_axes_list is None or bending_axis in bending_axes_list)
            and (support_list is None or support in support_list)
        )

        if not matches_filter:
            continue

        ##### Get section properties from AISC database ####
        r = (wide_flange_database[frame_name]['rx'] if bending_axis == 'x' else wide_flange_database[frame_name]['ry']) * inch

        ##### Reference interaction diagram ####
        analysis1_ALR_H = content[analysis1]['Original_ALR_H']
        analysis1_ALR_V = content[analysis1]['Original_ALR_V']

        analysis1_interaction_diagram = InteractionDiagram2d(
            analysis1_ALR_H,
            analysis1_ALR_V
        )

        #####################################################################
        # Compare analysis1 with every other analysis provided
        #####################################################################

        for analysis2 in analysis_types[1:]:

            analysis2_ALR_H = content[analysis2]['Original_ALR_H']
            analysis2_ALR_V = content[analysis2]['Original_ALR_V']

            theta_list = content[analysis2]['del2_over_del1_results']['Theta']
            del2_over_del1 = content[analysis2]['del2_over_del1_results']['del2_over_del1']
            del2_over_del1_capped=[min(x,cap) for x in del2_over_del1]

            analysis2_interaction_diagram = InteractionDiagram2d(
                analysis2_ALR_H,
                analysis2_ALR_V
            )

            errors = analysis1_interaction_diagram.compare_two(
                analysis2_interaction_diagram,
                theta_list,
                degrees=True
            )
            # print(theta_list)
            # print(del2_over_del1)
            # print(errors)
            # input("I am inside plot error vs del2 over del1")

            #################################################################
            # Store EVERY error with its corresponding del2/del1 value
            #################################################################

            for theta, error, ratio in zip(theta_list, errors, del2_over_del1_capped):

                results_data.append({
                    'frame': frame_id,
                    'section': frame_name,
                    'axis': bending_axis,
                    'stories': no_of_stories,
                    'height': storey_height,
                    'support': support,
                    'radius_of_gyration': r,
                    'comparison': f'{analysis1} vs {analysis2}',
                    'theta': theta,
                    'radial_error': error,
                    'del2_over_del1': ratio
                })

    #####################################################################
    # Create DataFrame
    #####################################################################

    results_df = pd.DataFrame(results_data)

    if results_df.empty:
        raise ValueError("No data matched the specified filters.")

    #####################################################################
    # Hover data
    #####################################################################

    hover_data_config = {
        'section': True,
        'axis': True,
        'stories': True,
        'height': ':.3f',
        'support': True,
        'radius_of_gyration': ':.3f',
        'comparison': True,
        'theta': ':.3f',
        'radial_error': ':.6f',
        'del2_over_del1': ':.6f'
    }

    #####################################################################
    # Create one scatter plot for each comparison and combine them
    #####################################################################

    fig = None

    colors = px.colors.qualitative.Plotly

    for i, analysis2 in enumerate(analysis_types[1:]):

        comparison_name = f'{analysis1} vs {analysis2}'

        comparison_df = results_df[
            results_df['comparison'] == comparison_name
        ]

        temp_fig = my_scatter_plot(
            comparison_df,
            x_col='del2_over_del1',
            y_col='radial_error',
            hover_data=hover_data_config,
            xlabel='Δ2 / Δ1',
            ylabel='Radial Error',
            title=f'Δ2/Δ1 vs Radial Error ({analysis1} vs Other Analyses)'
        )

        # Give each comparison its own color and legend name
        for trace in temp_fig.data:
            trace.name = comparison_name
            trace.showlegend = True
            trace.marker.color = colors[i % len(colors)]

        if fig is None:
            fig = temp_fig
        else:
            for trace in temp_fig.data:
                fig.add_trace(trace)

    #####################################################################
    # Save figure
    #####################################################################
    fig.add_hline(y=-0.05, line_dash="dot", line_color="black", line_width=1.5, opacity=0.9)
    os.makedirs(save_folder, exist_ok=True)

    base_filename = f"Error_vs_del2_over_del1_{analysis1}_vs_{'_'.join(analysis_types[1:])}"

    html_path = os.path.join(save_folder, f"{base_filename}.html")
    fig.write_html(html_path)
    print(f"Interactive plot saved to: {html_path}")

    png_path = os.path.join(save_folder, f"{base_filename}.png")
    fig.write_image(png_path, width=1000, height=600)
    print(f"Static PNG image saved to: {png_path}")

    fig.show()





def get_storey_height_list_for_a_section(
    section_name,
    bending_axis,
    slenderness_list,
    round_height=3
):


    if bending_axis == "x":
        r = wide_flange_database[section_name]["rx"] * inch
    elif bending_axis == "y":
        r = wide_flange_database[section_name]["ry"] * inch
    else:
        raise ValueError("bending_axis must be either 'x' or 'y'.")

    slenderness_list = np.array(slenderness_list, dtype=float)

    height_list = slenderness_list * r 

    return np.round(height_list, round_height).tolist()


def save_deformed_shape_frame(i, temporary_folder):
    temporary_folder = Path(temporary_folder)
    temporary_folder.mkdir(parents=True, exist_ok=True)

    # Create a separate figure only for this saved frame
    fig, ax = plt.subplots()

    opsvis.plot_defo(
        sfac=10,
        ax=ax,
        fmt_defo={
            "color": "blue",
            "linestyle": "solid",
            "linewidth": 3,
            "marker": "",
            "markersize": 1,
        },
        fmt_undefo={
            "color": "gray",
            "linestyle": "dashed",
            "linewidth": 2,
            "marker": "",
            "markersize": 1,
        },
    )

    ax.set_title(f"Step {i}")
    ax.axis("equal")

    frame_path = temporary_folder / f"{i:05d}.png"

    fig.savefig(
        frame_path,
        dpi=100,
        bbox_inches="tight",
    )

    # Close only the deformed-shape figure
    plt.close(fig)


def create_video_from_frames_and_clear_folder(
    temporary_folder,
    output_video_path,
    fps=10,
    pause_last_seconds=2
):

    temporary_folder = Path(temporary_folder)
    output_video_path = Path(output_video_path)

    # Make sure output folder exists
    output_video_path.parent.mkdir(parents=True, exist_ok=True)

    # Get and sort all png frames
    frame_files = sorted(temporary_folder.glob("*.png"))

    if len(frame_files) == 0:
        print("No frame images found. Video was not created.")
        return

    # Number of extra times to repeat the last frame
    last_frame_repeats = int(fps * pause_last_seconds)

    # Create video
    with imageio.get_writer(output_video_path, fps=fps) as writer:
        for frame_file in frame_files:
            image = imageio.imread(frame_file)
            writer.append_data(image)

        # Pause at final frame
        last_image = imageio.imread(frame_files[-1])

        for _ in range(last_frame_repeats):
            writer.append_data(last_image)

    print(f"Video saved at: {output_video_path}")

    # Delete all contents inside the temporary folder
    for item in temporary_folder.iterdir():
        if item.is_file() or item.is_symlink():
            item.unlink()
        elif item.is_dir():
            shutil.rmtree(item)

    print(f"All contents inside {temporary_folder} have been deleted.")


def save_analysis_history_figures(
    parent_folder,
    results,
    frame_name,
    analysis_type,
    ALR_H,
    ALR_V,
    step_size=None,
    tolerance=None,
    iterations=None,
    control_dir=None,
):
    """
    Plot six analysis-history curves in one 2-by-3 figure and save it.

    Parameters
    ----------
    parent_folder : str or Path
        Folder where the figure will be saved.

    results : object
        Analysis results containing the required response lists.

    frame_name : str
        Frame identifier.

    analysis_type : str
        Analysis type, such as GMNIA, GNA, or WC.

    ALR_H, ALR_V : float
        Horizontal and vertical axial load ratios.

    tolerance : float, optional
        Convergence tolerance used in the analysis.

    iterations : int, optional
        Maximum number of convergence iterations.

    control_dir : str, optional
        Control direction, such as "L" or "V".
    """

    parent_folder = Path(parent_folder)
    parent_folder.mkdir(
        parents=True,
        exist_ok=True,
    )

    fig, axes = plt.subplots(
        nrows=2,
        ncols=3,
        figsize=(16, 9),
    )

    axes = axes.ravel()

    plot_options = {
        "linewidth": 2,
        "marker": "o",
        "markersize": 3,
        "markevery": 1,
    }

    # 1. Displacement vs load ratio
    axes[0].plot(
        results.control_node_displacement,
        results.load_ratio,
        **plot_options,
    )
    axes[0].set_xlabel("Displacement at Control Node")
    axes[0].set_ylabel("Load Ratio λ")
    axes[0].set_title("Load Ratio vs Displacement")

    # 2. Load ratio vs eigenvalue
    axes[1].plot(
        results.load_ratio,
        results.lowest_eigenvalue,
        **plot_options,
    )
    axes[1].set_xlabel("Load Ratio λ")
    axes[1].set_ylabel("Lowest Eigenvalue")
    axes[1].set_title("Eigenvalue vs Load Ratio")

    # 3. Vertical reaction vs load ratio
    axes[2].plot(
        results.vertical_reaction,
        results.load_ratio,
        **plot_options,
    )
    axes[2].set_xlabel("Vertical Reaction")
    axes[2].set_ylabel("Load Ratio λ")
    axes[2].set_title("Load Ratio vs Vertical Reaction")

    # 4. Base shear vs load ratio
    axes[3].plot(
        results.base_shear,
        results.load_ratio,
        **plot_options,
    )
    axes[3].set_xlabel("Base Shear")
    axes[3].set_ylabel("Load Ratio λ")
    axes[3].set_title("Load Ratio vs Base Shear")

    # 5. Maximum strain vs load ratio
    axes[4].plot(
        results.absolute_maximum_strain,
        results.load_ratio,
        **plot_options,
    )
    axes[4].set_xlabel("Absolute Maximum Strain")
    axes[4].set_ylabel("Load Ratio λ")
    axes[4].set_title("Load Ratio vs Maximum Strain")

    # 6. Load ratio vs P-M-M interaction
    axes[5].plot(
        results.load_ratio,
        results.max_P_M_M_interaction,
        **plot_options,
    )
    axes[5].set_xlabel("Load Ratio λ")
    axes[5].set_ylabel("Maximum P-M-M Interaction")
    axes[5].set_title("P-M-M Interaction vs Load Ratio")

    # Common formatting
    for ax in axes:
        ax.grid(
            True,
            linestyle="--",
            linewidth=0.7,
            alpha=0.5,
        )

        ax.tick_params(
            axis="both",
            direction="in",
            labelsize=9,
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.margins(
            x=0.05,
            y=0.08,
        )

    # Main title
    title = f"{frame_name} — {analysis_type}"

    if control_dir is not None:
        title += f" — Control Direction: {control_dir}"

    fig.suptitle(
        title,
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )

    # Build information box
    information = [
        f"ALR_H: {ALR_H}",
        f"ALR_V: {ALR_V}",
    ]

    if tolerance is not None:
        information.append(
            f"Tolerance: {tolerance:.2e}"
        )

    if iterations is not None:
        information.append(
            f"Iterations: {iterations}"
        )

    if step_size is not None:
        information.append(
            f"No_of_steps: {step_size}"
        )

    information_text = "\n".join(information)

    fig.text(
        0.985,
        0.975,
        information_text,
        ha="right",
        va="top",
        fontsize=9,
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "white",
            "edgecolor": "black",
            "linewidth": 0.8,
            "alpha": 0.9,
        },
    )

    fig.tight_layout(
        rect=[0, 0, 1, 0.91],
        h_pad=2,
        w_pad=2,
    )

    # Build filename dynamically
    filename_parts = [
        f"ALR_H_{ALR_H}",
        f"ALR_V_{ALR_V}",
        # str(frame_name),
        # str(analysis_type),
    ]

    if control_dir is not None:
        filename_parts.append(
            f"control_{control_dir}"
        )

    if tolerance is not None:
        tolerance_text = (
            f"{tolerance:.2e}"
            .replace("+", "")
        )

        filename_parts.append(
            f"tol_{tolerance_text}"
        )

    if iterations is not None:
        filename_parts.append(
            f"iter_{iterations}"
        )

    if step_size is not None:
        filename_parts.append(
            f"No_of_steps_{step_size}"
        )

    filename = parent_folder / (
        "_".join(filename_parts) + ".png"
    )

    fig.savefig(
        filename,
        dpi=300,
        bbox_inches="tight",
    )

    # Close only this figure
    plt.close(fig)

def check_load_ratio_problem(
    results,
    frame_name,
    analysis_type,
    ALR_V,
    ALR_H,
    txt_file_path="problematic_cases.txt",
    window=11,
    robust_z_limit=6.0,
    minimum_consecutive_flags=1,
):
    """
    Detect abrupt changes in the slope of the load-displacement curve.

    The response slope is:

        slope = delta(load ratio) / delta(displacement)

    An abrupt change is identified by comparing the change in slope
    against the recent local distribution using median absolute
    deviation (MAD).
    """

    load_ratio = np.asarray(
        results.load_ratio,
        dtype=float,
    )

    displacement = np.asarray(
        results.control_node_displacement,
        dtype=float,
    )

    if len(load_ratio) != len(displacement):
        raise ValueError(
            "load_ratio and control_node_displacement "
            "must have equal lengths."
        )

    if len(load_ratio) < max(window + 2, 6):
        return "Not Problematic"

    window = max(int(window), 5)

    delta_lambda = np.diff(load_ratio)
    delta_u = np.diff(displacement)

    # Avoid division by zero or almost-zero displacement increments
    displacement_scale = max(
        np.max(np.abs(displacement)),
        1.0,
    )

    displacement_tolerance = (
        100
        * np.finfo(float).eps
        * displacement_scale
    )

    valid_increment = (
        np.abs(delta_u) > displacement_tolerance
    )

    slopes = np.full(
        len(delta_lambda),
        np.nan,
        dtype=float,
    )

    slopes[valid_increment] = (
        delta_lambda[valid_increment]
        / delta_u[valid_increment]
    )

    # Change in slope: a discrete measure of curvature
    slope_changes = np.diff(slopes)

    flagged_slope_change_indices = []
    details = []

    for i in range(window, len(slope_changes)):
        current_change = slope_changes[i]

        if not np.isfinite(current_change):
            continue

        # Use only earlier values so the detector can also work online
        local_changes = slope_changes[i - window:i]

        local_changes = local_changes[
            np.isfinite(local_changes)
        ]

        if len(local_changes) < 4:
            continue

        local_median = np.median(local_changes)

        mad = np.median(
            np.abs(local_changes - local_median)
        )

        robust_scale = 1.4826 * mad

        # Fallback when the nearby response is nearly perfectly smooth
        if robust_scale == 0:
            q1, q3 = np.percentile(
                local_changes,
                [25, 75],
            )

            robust_scale = (q3 - q1) / 1.349

        numerical_floor = (
            1e-10
            * max(
                np.median(np.abs(slopes[np.isfinite(slopes)])),
                1.0,
            )
        )

        robust_scale = max(
            robust_scale,
            numerical_floor,
        )

        robust_z_score = (
            abs(current_change - local_median)
            / robust_scale
        )

        if robust_z_score > robust_z_limit:
            # slope_changes[i] represents the change between
            # slopes[i] and slopes[i + 1].
            #
            # slopes[i + 1] represents the result increment ending
            # at result index i + 2.
            result_index = i + 2

            flagged_slope_change_indices.append(
                result_index
            )

            details.append(
                {
                    "index": result_index,
                    "previous_load_ratio": load_ratio[
                        result_index - 1
                    ],
                    "current_load_ratio": load_ratio[
                        result_index
                    ],
                    "previous_displacement": displacement[
                        result_index - 1
                    ],
                    "current_displacement": displacement[
                        result_index
                    ],
                    "previous_slope": slopes[i],
                    "current_slope": slopes[i + 1],
                    "slope_change": current_change,
                    "expected_slope_change": local_median,
                    "robust_z_score": robust_z_score,
                }
            )

    # Optionally require multiple adjacent flags
    if minimum_consecutive_flags > 1:
        retained_indices = []

        for j, index in enumerate(
            flagged_slope_change_indices
        ):
            nearby_count = sum(
                abs(index - other_index)
                <= minimum_consecutive_flags
                for other_index
                in flagged_slope_change_indices
            )

            if nearby_count >= minimum_consecutive_flags:
                retained_indices.append(index)

        retained_set = set(retained_indices)

        details = [
            detail
            for detail in details
            if detail["index"] in retained_set
        ]

        flagged_slope_change_indices = retained_indices

    status = (
        "Problematic"
        if flagged_slope_change_indices
        else "Not Problematic"
    )

    if status == "Problematic":
        txt_file_path = Path(txt_file_path)
        txt_file_path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        with open(
            txt_file_path,
            "a",
            encoding="utf-8",
        ) as file:
            file.write(
                f"\nFrame: {frame_name}, "
                f"Analysis: {analysis_type}, "
                f"ALR_V: {ALR_V}, "
                f"ALR_H: {ALR_H}\n"
            )

            for detail in details:
                file.write(
                    f"    Index: {detail['index']}, "
                    f"Load ratio: "
                    f"{detail['previous_load_ratio']:.8g} -> "
                    f"{detail['current_load_ratio']:.8g}, "
                    f"Displacement: "
                    f"{detail['previous_displacement']:.8g} -> "
                    f"{detail['current_displacement']:.8g}, "
                    f"Slope: "
                    f"{detail['previous_slope']:.8g} -> "
                    f"{detail['current_slope']:.8g}, "
                    f"Slope change: "
                    f"{detail['slope_change']:.8g}, "
                    f"Robust score: "
                    f"{detail['robust_z_score']:.3f}\n"
                )

