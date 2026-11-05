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
from filelock import FileLock
import runpy
from matplotlib.lines import Line2D
from tqdm import tqdm
from scipy.interpolate import griddata

CONFIG_PATH = Path(__file__).resolve().parent / "Columns_config.py"
CONFIG_LOCK = FileLock(str(CONFIG_PATH) + ".lock", timeout=120)

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

    file_path = CONFIG_PATH

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
        raise ValueError("Frame_Info dictionary was not found in Columns_config.py.")
    content_bytes = content.encode("utf-8")
    lines = content.splitlines(keepends=True)

    start_index = (
        sum(len(line.encode("utf-8")) for line in lines[:frame_info_node.lineno - 1])
        + frame_info_node.col_offset
    )

    end_index = (
        sum(len(line.encode("utf-8")) for line in lines[:frame_info_node.end_lineno - 1])
        + frame_info_node.end_col_offset
    )
    updated_frame_info = pformat(Frame_Info, sort_dicts=False, width=120).encode("utf-8")

    updated_content = (
        content_bytes[:start_index]
        + updated_frame_info
        + content_bytes[end_index:]
    )

    compile(updated_content, str(file_path), "exec")

    temp_path = None

    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=file_path.parent,
            suffix=".tmp",
            delete=False
        ) as file:

            temp_path = Path(file.name)
            file.write(updated_content)
            file.flush()
            os.fsync(file.fileno())

        os.replace(temp_path, file_path)

    finally:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink()

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
        "Base_Wind_load": 1000 * KN,
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


    with CONFIG_LOCK:

        # Reload the latest dictionary from Columns_config.py
        latest_frame_info = runpy.run_path(str(CONFIG_PATH))["Frame_Info"]

        # Update the existing dictionary in memory
        Frame_Info.clear()
        Frame_Info.update(latest_frame_info)

        if code not in Frame_Info:
            Frame_Info[code] = config
            _save_frame_info_to_config()
            print(f"Created new configuration: {code}", flush=True)

        elif Frame_Info[code] != config:
            Frame_Info[code] = config
            _save_frame_info_to_config()
            print(f"Updated existing configuration: {code}", flush=True)

        else:
            print(f"Configuration already exists and is unchanged: {code}", flush=True)

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
    analysis_types=None
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
        analyses = [key for key, value in content.items() 
                    if key != "data" and isinstance(value, dict) 
                    and (analysis_types is None or key in analysis_types)]
        
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

            filename = f"interaction_diagrams.png"
            save_path = os.path.join(save_folder, filename)

            fig.savefig(save_path, dpi=300, bbox_inches="tight")

            print(f"[{file_number}/{len(json_files)}] Saved to: {save_path}")

        if show:
            plt.show()
        else:
            plt.close(fig)

def plot_interaction_diagrams_for_all_analyses_with_load_combo_info(json_file_path, use_original=True, save=True, show=True, analysis_types=None, load_combo_analysis=None, load_combo_L_over_D_ratio=None, load_combo_E_over_D_ratio=None):

    Folder = Path(json_file_path)

    if not Folder.exists():
        raise FileNotFoundError(f"The specified path does not exist: {Folder}")

    if not Folder.is_dir():
        raise NotADirectoryError(f"The specified path is not a directory: {Folder}")

    json_files = sorted(file_path for file_path in Folder.iterdir() if file_path.is_file() and file_path.suffix.lower() == ".json")

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

        analyses = [key for key, value in content.items() if key != "data" and isinstance(value, dict) and (analysis_types is None or key in analysis_types)]

        if not analyses:
            print(f"Skipping {file_path.name}: no analyses were found.")
            continue

        fig, ax = plt.subplots(figsize=(8, 6))
        plotted_analyses = []
        line_styles = ["-", "--", "-.", ":"]

        for analysis_index, analysis in enumerate(analyses):

            analysis_data = content[analysis]

            if alr_h_key not in analysis_data or alr_v_key not in analysis_data:
                print(f"Skipping analysis '{analysis}' in {file_path.name}: missing {alr_h_key} or {alr_v_key}.")
                continue

            ALR_H = analysis_data[alr_h_key]
            ALR_V = analysis_data[alr_v_key]

            if ALR_H is None or ALR_V is None:
                print(f"Skipping analysis '{analysis}' in {file_path.name}: interaction data are None.")
                continue

            if len(ALR_H) == 0 or len(ALR_V) == 0:
                print(f"Skipping analysis '{analysis}' in {file_path.name}: interaction data are empty.")
                continue

            if len(ALR_H) != len(ALR_V):
                print(f"Skipping analysis '{analysis}' in {file_path.name}: ALR_H and ALR_V have different lengths.")
                continue

            interaction_diagram = InteractionDiagram2d(ALR_H, ALR_V)

            plt.sca(ax)
            interaction_diagram.plot(label=analysis, linestyle=line_styles[analysis_index % len(line_styles)], linewidth=1.8, alpha=0.9)

            plotted_analyses.append(analysis)

        if not plotted_analyses:
            print(f"Skipping {file_path.name}: no valid interaction curves were found.")
            plt.close(fig)
            continue

        #### Plot load combination points ####

        if load_combo_analysis is not None:

            if load_combo_analysis not in content:
                print(f"Load combination analysis '{load_combo_analysis}' was not found in {file_path.name}.")

            else:
                load_combo_cases = content[load_combo_analysis].get("load_combo_details", {})
                load_combo_points_plotted = False

                for ratio_key, load_combo_case in load_combo_cases.items():
                    stored_L_over_D_ratio = load_combo_case.get("L_over_D_ratio")
                    stored_E_over_D_ratio = load_combo_case.get("E_over_D_ratio")

                    matches_L_over_D = load_combo_L_over_D_ratio is None or stored_L_over_D_ratio == load_combo_L_over_D_ratio
                    matches_E_over_D = load_combo_E_over_D_ratio is None or stored_E_over_D_ratio == load_combo_E_over_D_ratio

                    if not (matches_L_over_D and matches_E_over_D):
                        continue

                    load_combo_results = load_combo_case["load_combo_results"]

                    for load_combo_name, load_combo_result in load_combo_results.items():
                        if load_combo_name in ["Dmax","governing_load_combo"]:
                            continue
                        horizontal_load = load_combo_result["horizontal_load"]
                        vertical_load = load_combo_result["vertical_load"]
                        governing = load_combo_result["governing"]

                        if governing:
                            ax.scatter(horizontal_load, vertical_load, s=50, marker="*", color="red", edgecolor="black", linewidth=0.8, zorder=10, label="_nolegend_")
                        else:
                            ax.scatter(horizontal_load, vertical_load, s=10, marker="o", color="blue", zorder=9, label="_nolegend_")

                        load_combo_points_plotted = True

                if load_combo_points_plotted:
                    ax.scatter([], [], s=140, marker="*", color="red", edgecolor="black", linewidth=0.8, label="Controlling load combination")
        ax.set_xlabel("ALR_H")
        ax.set_ylabel("ALR_V")
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        ax.legend(title="Analysis / Load Combination", loc="best")

        if no_of_stories == 1:
            story_label = "1 story"
        elif isinstance(no_of_stories, (int, float)):
            story_label = f"{no_of_stories} stories"
        else:
            story_label = str(no_of_stories)

        title_lines = [f"{section_name}, axis = {bending_axis}", f"{story_label}, h = {storey_height} m"]

        if material is not None:
            title_lines.append(f"Material = {material}")

        if load_combo_analysis is not None:
            L_D_label = "all" if load_combo_L_over_D_ratio is None else str(load_combo_L_over_D_ratio)
            E_D_label = "all" if load_combo_E_over_D_ratio is None else str(load_combo_E_over_D_ratio)
            title_lines.append(f"Load combinations: {load_combo_analysis}, L/D = {L_D_label}, E/D = {E_D_label}")

        ax.set_title("\n".join(title_lines))
        fig.tight_layout()

        save_path = None

        if save:
            save_folder = os.path.join("Column_Results/Analysis_History", f"{frame_id}")
            os.makedirs(save_folder, exist_ok=True)
            L_D_label = "all" if load_combo_L_over_D_ratio is None else str(load_combo_L_over_D_ratio)
            E_D_label = "all" if load_combo_E_over_D_ratio is None else str(load_combo_E_over_D_ratio)
            filename = f"interaction_diagrams_L_D_{L_D_label}_E_D_{E_D_label}.png"
            save_path = os.path.join(save_folder, filename)
            fig.savefig(save_path, dpi=300, bbox_inches="tight")
            print(f"[{file_number}/{len(json_files)}] Saved to: {save_path}")

        if show:
            plt.show()
        else:
            plt.close(fig)


def plot_governing_D_max_vs_E_over_D_for_all_analyses(json_file_path, save=True, show=True, analysis_types=None, L_over_D_ratios=None, E_over_D_ratios=None):

    load_combo_colors = {"1.4D": "red", "1.2D+1.6L": "blue", "1.2D+L+E": "green", "0.9D+E": "orange"}
    linestyles = ['-', '--', '-.', ':']

    for file_path in sorted(Path(json_file_path).iterdir()):
        if not file_path.is_file() or file_path.suffix.lower() != ".json": continue
        with open(file_path, "r", encoding="utf-8") as file:
            content = json.load(file)

        data = content.get("data", {})
        frame_id = data.get("frame_id", file_path.stem)
        section_name = data.get("column_section_name", "Unknown section")
        bending_axis = data.get("bending_axes", data.get("bending_axis", "Unknown axis"))
        storey_height = data.get("storey_height", "Unknown height")
        no_of_stories = data.get("no_of_stories", "Unknown number of stories")
        material = data.get("Material", data.get("Steel_Grade", None))

        analyses = [key for key, value in content.items() if key != "data" and isinstance(value, dict) and "load_combo_details" in value and (analysis_types is None or key in analysis_types)]

        if not analyses:
            print(f"Skipping {file_path.name}: no analyses with load combination details were found.")
            continue

        save_folder = Path("Column_Results/Analysis_History") / str(frame_id) / "D_max_vs_E_over_D_plots"
        if save: save_folder.mkdir(parents=True, exist_ok=True)

        for L_over_D in L_over_D_ratios:
            RESULTS = {}
            RESULTS['L_over_D'] = L_over_D
            print(f"Analysis for {L_over_D} running")
            fig, ax = plt.subplots(figsize=(8, 6))
            analysis_handles = []
            for i, analysis in enumerate(analyses):
                results = get_Dmax_EoverD_LoverD(content, analysis_type=analysis, L_over_D_ratios=[L_over_D], E_over_D_ratios=E_over_D_ratios)
                RESULTS[analysis] = {}
                RESULTS[analysis]['E_over_D'] = [v['E_over_D_ratio'] for v in results.values()]
                RESULTS[analysis]['Dmax'] = [v['Dmax'] for v in results.values()]
                RESULTS[analysis]['Governing_load_combo'] = [v['Governing_Load_Combo'] for v in results.values()]

                E_over_D = np.array(RESULTS[analysis]['E_over_D'])
                Dmax = np.array(RESULTS[analysis]['Dmax'])
                governing_combos = RESULTS[analysis]['Governing_load_combo']
                if len(E_over_D) == 0: continue
                sorted_indices = np.argsort(E_over_D)
                line, = ax.plot(E_over_D[sorted_indices], Dmax[sorted_indices], linestyle=linestyles[i % len(linestyles)], linewidth=1.8, label=analysis)
                ax.scatter(E_over_D, Dmax, c=[load_combo_colors.get(combo, 'gray') for combo in governing_combos], s=15, zorder=3)

                analysis_handles.append(line)

            ax.set_xlabel('E/D Ratio')
            ax.set_ylabel('D_max')
            ax.set_title(f'L/D = {L_over_D}')
            ax.grid(True, alpha=0.3)

            combo_handles = [Line2D([0], [0], marker='o', color='none', markerfacecolor=color, markersize=8, linestyle='None', label=combo) for combo, color in load_combo_colors.items()]
            if analysis_handles:
                legend1 = ax.legend(handles=analysis_handles, loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=len(analysis_handles), frameon=False)
                ax.add_artist(legend1)

            ax.legend(handles=combo_handles, loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=len(combo_handles), frameon=False)
            fig.subplots_adjust(bottom=0.30)

            if save: fig.savefig(save_folder / f'Dmax_vs_EoverD_LoverD_{L_over_D:.4f}.png', dpi=300, bbox_inches='tight')
            if show: plt.show()
            plt.close(fig)


def plot_theta_vs_del2_over_del1_for_all_analyses(json_file_path, cap=3, save=True, show=True):

    Folder = Path(json_file_path)
    results = []
    line_styles = ['-', '--', '-.', ':']

    for file_path in tqdm(sorted(Folder.glob("*.json")), desc="Processing Frames", unit="frame"):

        with open(file_path, "r", encoding="utf-8") as file:
            content = json.load(file)

        # Frame details
        data = content['data']
        frame_id = data['Frame_id']
        section_name = data['column_section']['common_and_exceptions']['common'][0]
        bending_axis = data['column_section']['common_and_exceptions']['common'][1]
        storey_height = data['storey_height'][0]
        no_of_stories = len(data['storey_height'])

        plt.figure(figsize=(8, 6))
        plotted_analyses = []

        # Plot each analysis
        analyses = [key for key, value in content.items() if key != 'data' and isinstance(value, dict)]

        for i, analysis in enumerate(analyses):

            if 'del2_over_del1_results' not in content[analysis]:
                continue

            theta = content[analysis]['del2_over_del1_results']['Theta']
            del2_over_del1 = content[analysis]['del2_over_del1_results']['del2_over_del1']
            del2_over_del1_capped = [min(x, cap) for x in del2_over_del1]

            plt.plot(theta, del2_over_del1_capped, label=analysis, linestyle=line_styles[i % len(line_styles)], linewidth=1.5)
            plotted_analyses.append(analysis)

        if not plotted_analyses:
            plt.close()
            continue

        # Plot formatting
        plt.xlabel(r'$\theta$ (degrees)')
        plt.ylabel(r'$\Delta_2/\Delta_1$')
        plt.title(f'{section_name}, axis = {bending_axis}\n{no_of_stories} stories, h = {storey_height} m')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        # Save figure
        save_path = None

        if save:
            save_folder = os.path.join("Column_Results", "Analysis_History", str(frame_id))
            os.makedirs(save_folder, exist_ok=True)
            save_path = os.path.join(save_folder, f"{frame_id}_theta_vs_del2_over_del1.png")
            plt.savefig(save_path, dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        else:
            plt.close()

        results.append({'frame_id': frame_id, 'analyses': plotted_analyses, 'save_path': save_path})

    return results

def plot_theta_vs_error_all_analyses(json_file_path, benchmark_analysis, test_analyses):

    Folder = Path(json_file_path)
    theta_list = np.linspace(0, 90, 181)

    for file_path in tqdm(sorted(Folder.glob("*.json")), desc="Processing Frames", unit="frame"):

        with open(file_path, "r", encoding="utf-8") as file:
            content = json.load(file)

        frame_id = content['data']['Frame_id']

        # Benchmark interaction diagram
        benchmark_ALR_H = content[benchmark_analysis]['Original_ALR_H']
        benchmark_ALR_V = content[benchmark_analysis]['Original_ALR_V']
        benchmark_interaction_diagram = InteractionDiagram2d(benchmark_ALR_H, benchmark_ALR_V)

        plt.figure(figsize=(8, 6))

        # Compare benchmark with each test analysis
        for test_analysis in test_analyses:

            test_ALR_H = content[test_analysis]['Original_ALR_H']
            test_ALR_V = content[test_analysis]['Original_ALR_V']
            test_interaction_diagram = InteractionDiagram2d(test_ALR_H, test_ALR_V)

            errors = benchmark_interaction_diagram.compare_two(test_interaction_diagram, theta_list, degrees=True)

            plt.plot(theta_list, np.array(errors) * 100, label=test_analysis, linewidth=1.5)

        # Plot formatting
        plt.axhline(0, color='black', linewidth=1)
        plt.xlabel(r'$\theta$ (degrees)')
        plt.ylabel('Radial Error (%)')
        plt.title(f'{frame_id}: {benchmark_analysis} vs Test Analyses')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        save_folder = os.path.join("Column_Results", "Analysis_History", str(frame_id))
        os.makedirs(save_folder, exist_ok=True)
        plt.savefig(os.path.join(save_folder, f"{frame_id}_theta_vs_error.png"), dpi=300, bbox_inches='tight')
        plt.close()


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
            size=6,
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
        frame_id = content['data']['Frame_id']
        frame_name = content['data']['column_section']['common_and_exceptions']['common'][0]  ### !!!! This works only for MFColumn Results
        bending_axis = content['data']['column_section']['common_and_exceptions']['common'][1]### !!!! This works only for MFColumn Results
        storey_height = content['data']["storey_height"][0]
        no_of_stories = len(content['data']["storey_height"])
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

def plot_slenderness_ratio_vs_del2_over_del1_vs_error_contour_plots(
    json_file_path, benchmark_analysis, test_analysis,
    frame_names_list=None, no_of_stories_list=None,
    bending_axes_list=None, support_list=None,
    min_error=-6, max_error=-1, save=True, show=False):

    save_folder = os.path.join("Column_Results", "slenderness_ratio_vs_del2_over_del1_vs_error_contour_plots")

    Folder = Path(json_file_path)
    slenderness_ratio = []
    del2_over_del1_list = []
    error_list = []

    json_files = list(Folder.glob("*.json"))

    for file_path in tqdm(json_files, desc="Processing Frames", unit="frame"):

        tqdm.write(f"Parsing: {file_path.name}")

        with open(file_path, "r", encoding="utf-8") as file:
            content = json.load(file)

        # Read structural parameters
        frame_name = content['data']['column_section']['common_and_exceptions']['common'][0]
        bending_axis = content['data']['column_section']['common_and_exceptions']['common'][1]
        storey_height = content['data']['storey_height'][0]
        no_of_stories = len(content['data']['storey_height'])
        support = content['data']['support']

        # Apply filters
        matches_filter = (
            (frame_names_list is None or frame_name in frame_names_list)
            and (no_of_stories_list is None or no_of_stories in no_of_stories_list)
            and (bending_axes_list is None or bending_axis in bending_axes_list)
            and (support_list is None or support in support_list)
        )

        if not matches_filter:
            continue

        # Calculate slenderness ratio
        r = (wide_flange_database[frame_name]['rx'] if bending_axis == 'x' else wide_flange_database[frame_name]['ry']) * inch
        slenderness_ratio.append(storey_height / r)

        # Benchmark interaction diagram
        benchmark_ALR_H = content[benchmark_analysis]['Original_ALR_H']
        benchmark_ALR_V = content[benchmark_analysis]['Original_ALR_V']
        benchmark_interaction_diagram = InteractionDiagram2d(benchmark_ALR_H, benchmark_ALR_V)

        # Test interaction diagram
        test_ALR_H = content[test_analysis]['Original_ALR_H']
        test_ALR_V = content[test_analysis]['Original_ALR_V']
        test_interaction_diagram = InteractionDiagram2d(test_ALR_H, test_ALR_V)

        # Get theta and displacement amplification
        test_results = content[test_analysis]['del2_over_del1_results']
        angles = test_results['Theta']
        del2_over_del1 = test_results['del2_over_del1']

        # Calculate radial errors
        errors = benchmark_interaction_diagram.compare_two(test_interaction_diagram, angles, degrees=True)

        del2_over_del1_list.append(del2_over_del1)
        error_list.append(errors)

    # Flatten data
    if not slenderness_ratio:
        print("No frames matched the specified filters.")
        return

    X = np.concatenate([np.full(len(values), slenderness) for slenderness, values in zip(slenderness_ratio, del2_over_del1_list)])
    Y = np.concatenate(del2_over_del1_list)
    Z = np.concatenate(error_list) * 100

    # Filter based on radial error only
    if min_error >= max_error:
        raise ValueError("min_error must be smaller than max_error.")

    valid = np.isfinite(X) & np.isfinite(Y) & np.isfinite(Z) & (Z >= min_error) & (Z <= max_error)
    X, Y, Z = X[valid], Y[valid], Z[valid]

    if len(X) < 3 or len(np.unique(X)) < 2 or len(np.unique(Y)) < 2:
        print("Insufficient data for contour plotting after filtering.")
        return

    # Create interpolation grid
    x_grid = np.linspace(np.min(X), np.max(X), 200)
    y_grid = np.linspace(np.min(Y), np.max(Y), 200)
    Xg, Yg = np.meshgrid(x_grid, y_grid)
    Zg = griddata((X, Y), Z, (Xg, Yg), method='linear')

    # Contour plot
    plt.figure(figsize=(8, 6))
    levels = np.linspace(min_error, max_error, 21)
    cs = plt.contour(Xg, Yg, Zg, levels=levels, cmap='viridis')
    plt.clabel(cs, inline=True, fontsize=10, fmt='%g')

    plt.xlabel('Member Slenderness, L/r')
    plt.ylabel(r'Ratio of Second-Order Drift to First-Order Drift, $\Delta_2/\Delta_1$')
    plt.title(f'{benchmark_analysis} vs {test_analysis}')
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(np.min(Y), np.max(Y))
    plt.tight_layout()

    # Save figure
    if save:
        os.makedirs(save_folder, exist_ok=True)
        plt.savefig(os.path.join(save_folder, f'{benchmark_analysis}_vs_{test_analysis}_contour_lines.png'), dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close()

def get_governing_load_combo_details(json_file_path, analysis_type, L_over_D_ratios, E_over_D_ratios,
                                     frame_names_list=None, storey_heights_list=None, no_of_stories_list=None, bending_axes_list=None, support_list=None):

    for file_path in Path(json_file_path).iterdir():

        # Find and load json files
        if not file_path.is_file() or file_path.suffix.lower() != ".json":
            continue

        with open(file_path, "r", encoding="utf-8") as file:
            content = json.load(file)

        # Filters to run only selected frames
        frame_id = content["data"]["Frame_id"]
        frame_name = content["data"]["column_section_name"]
        bending_axis = content["data"]["bending_axes"]
        storey_height = content["data"]["storey_height"][0]
        no_of_stories = content["data"]["no_of_stories"]
        support = content["data"]["support"]

        matches_filter = ((frame_names_list is None or frame_name in frame_names_list) 
                    and (storey_heights_list is None or storey_height in storey_heights_list)
                    and (no_of_stories_list is None or no_of_stories in no_of_stories_list) 
                    and (bending_axes_list is None or bending_axis in bending_axes_list) 
                    and (support_list is None or support in support_list))

        if not matches_filter:
            continue

        # Reteive the interaction diagram
        analysis1_ALR_H = content[analysis_type]["Original_ALR_H"]
        analysis1_ALR_V = content[analysis_type]["Original_ALR_V"]
        analysis1_interaction_diagram = InteractionDiagram2d(analysis1_ALR_H, analysis1_ALR_V)

        # Initilize results dictionary
        content[analysis_type]["load_combo_details"] = {}

        # Computes results
        for L_over_D_ratio in L_over_D_ratios:
            for E_over_D_ratio in E_over_D_ratios:
                load_combo_results = find_the_intersection_of_load_combo_with_interaction_curves(analysis1_interaction_diagram, L_over_D_ratio, E_over_D_ratio)
                load_combo_details = {"L_over_D_ratio": L_over_D_ratio, "E_over_D_ratio": E_over_D_ratio, "load_combo_results": load_combo_results}
                ratio_key = f"L_D_{L_over_D_ratio}__E_D_{E_over_D_ratio}"
                content[analysis_type]["load_combo_details"][ratio_key] = load_combo_details
                print(f"Prepared: {file_path.name} | {ratio_key}")

        # Save json file
        with open(file_path, "w", encoding="utf-8") as file:
            json.dump(content, file, indent=4, default=lambda x: int(x) if isinstance(x, np.integer) else x)

        print(f"Updated file: {file_path.name}")



def get_Dmax_EoverD_LoverD(json_file_content, analysis_type, L_over_D_ratios, E_over_D_ratios):


    # Retrieve the interaction diagram
    analysis1_ALR_H = json_file_content[analysis_type]["Original_ALR_H"]
    analysis1_ALR_V = json_file_content[analysis_type]["Original_ALR_V"]
    analysis1_interaction_diagram = InteractionDiagram2d(analysis1_ALR_H, analysis1_ALR_V)

    # Initialize results dictionary
    results = {}
    x = 0

    # Compute results
    for L_over_D_ratio in L_over_D_ratios:
        for E_over_D_ratio in E_over_D_ratios:
            load_combo_results = find_the_intersection_of_load_combo_with_interaction_curves(
                analysis1_interaction_diagram,
                L_over_D_ratio,
                E_over_D_ratio
            )
            results[x] = {
                'L_over_D_ratio': L_over_D_ratio,
                'E_over_D_ratio': E_over_D_ratio,
                'Dmax': load_combo_results['Dmax'],
                'Governing_Load_Combo': load_combo_results['governing_load_combo']
            }
            x += 1
    return results

def find_the_intersection_of_load_combo_with_interaction_curves(interaction_diagram, L_over_D_ratio, E_over_D_ratio):

    load_combinations = {
        "1.4D":      [0, 1.4],
        "1.2D+1.6L": [0, 1.2 + 1.6 * L_over_D_ratio],
        "1.2D+L+E":  [E_over_D_ratio, 1.2 + L_over_D_ratio],
        "0.9D+E":    [E_over_D_ratio, 0.9],
    }

    # Evaluate each load combination
    ratio_of_load_values_to_strengths = {}
    for load_combo_name, load_factors in load_combinations.items():
        horizontal_load, vertical_load = load_factors
        ratio = interaction_diagram.check_points(horizontal_load, vertical_load)
        ratio_of_load_values_to_strengths[load_combo_name] = ratio

    # Determine the governing load combination
    governing_load_combo = max(ratio_of_load_values_to_strengths, key=ratio_of_load_values_to_strengths.get)
    governing_ratio = ratio_of_load_values_to_strengths[governing_load_combo]

    # Save results
    results = {}
    results['Dmax'] = 1 / governing_ratio
    results['governing_load_combo'] = governing_load_combo   
    # for load_combo_name, load_factors in load_combinations.items():
    #     horizontal_load, vertical_load = load_factors  
    #     results[load_combo_name] = {"horizontal_load": horizontal_load/governing_ratio, "vertical_load": vertical_load/governing_ratio, 'governing': load_combo_name == governing_load_combo}
    
    return results

def plot_LoverD_EoverD_Dmaxerror_contour_plots(json_file_path,benchmark_analysis,test_analysis,L_over_D_ratios, E_over_D_ratios,
                                     frame_names_list=None, storey_heights_list=None, no_of_stories_list=None, bending_axes_list=None, support_list=None,save=True,show=False):

    json_files = [f for f in Path(json_file_path).iterdir() if f.is_file() and f.suffix.lower() == ".json"]

    for file_path in tqdm(json_files, desc="Processing Frames", unit="frame"):

        # Find and load json files
        if not file_path.is_file() or file_path.suffix.lower() != ".json":
            continue

        with open(file_path, "r", encoding="utf-8") as file:
            content = json.load(file)

        # Filters to run only selected frames
        frame_id = content['data']['Frame_id']
        frame_name = content['data']['column_section']['common_and_exceptions']['common'][0]  ### !!!! This works only for MFColumn Results
        bending_axis = content['data']['column_section']['common_and_exceptions']['common'][1]### !!!! This works only for MFColumn Results
        storey_height = content['data']["storey_height"][0]
        no_of_stories = len(content['data']["storey_height"])
        support = content['data']["support"]

        matches_filter = ((frame_names_list is None or frame_name in frame_names_list) 
                    and (storey_heights_list is None or storey_height in storey_heights_list)
                    and (no_of_stories_list is None or no_of_stories in no_of_stories_list) 
                    and (bending_axes_list is None or bending_axis in bending_axes_list) 
                    and (support_list is None or support in support_list))

        if not matches_filter:
            continue

        L_over_D_all = []
        E_over_D_all = []
        Dmax_error_all = []

        for L_over_D in tqdm(L_over_D_ratios, desc=f"Frame {frame_id}", unit="ratio", leave=False):
            results_benchmark=get_Dmax_EoverD_LoverD(content, benchmark_analysis, L_over_D_ratios=[L_over_D], E_over_D_ratios=E_over_D_ratios)
            results_test=get_Dmax_EoverD_LoverD(content, test_analysis, L_over_D_ratios=[L_over_D], E_over_D_ratios=E_over_D_ratios)
            E_over_D = [v['E_over_D_ratio'] for v in results_benchmark.values()]
            Dmax_benchmark = [v['Dmax'] for v in results_benchmark.values()]
            Dmax_test = [v['Dmax'] for v in results_test.values()]
            Dmax_error = [(benchmark - test) / benchmark * 100 for benchmark, test in zip(Dmax_benchmark, Dmax_test)]
            
            # Store data for contour plot
            L_over_D_all.extend([L_over_D] * len(E_over_D))
            E_over_D_all.extend(E_over_D)
            Dmax_error_all.extend(Dmax_error)

        E_grid, L_grid = np.meshgrid(E_over_D_ratios, L_over_D_ratios)
        Error_grid = np.array(Dmax_error_all).reshape(len(L_over_D_ratios), len(E_over_D_ratios))

        save_folder = Path("Column_Results") / "Analysis_History" / str(frame_id) / "Contour plots"
        if save:
            save_folder.mkdir(parents=True, exist_ok=True)

        fig, ax = plt.subplots(figsize=(8, 6))

        contours = ax.contour(E_grid, L_grid, Error_grid, levels=10, colors='black', linewidths=0.8)
        ax.clabel(contours, inline=True, fontsize=8, fmt='%.1f%%')

        if np.nanmin(Error_grid) <= -5 <= np.nanmax(Error_grid):
            contour_minus5 = ax.contour(E_grid, L_grid, Error_grid, levels=[-5], colors='red', linewidths=2)
            ax.clabel(contour_minus5, inline=True, fontsize=9, fmt='%.1f%%')

        ax.set_xlabel('E/D Ratio')
        ax.set_ylabel('L/D Ratio')
        ax.set_title(f'Dmax Error Contour: {benchmark_analysis} vs {test_analysis}')
        ax.grid(alpha=0.2)
        plt.tight_layout()

        if save:
            fig.savefig(save_folder / f'Dmax_Error_Contour_{benchmark_analysis}_vs_{test_analysis}.png', dpi=300, bbox_inches='tight')
        if show:
            plt.show()
        plt.close(fig)



#################################################################################
        fig, ax = plt.subplots(figsize=(8, 6))

        mesh = ax.pcolormesh(E_grid, L_grid, Error_grid, shading='auto', cmap='coolwarm')

        for i, L_over_D in enumerate(L_over_D_ratios):
            for j, E_over_D in enumerate(E_over_D_ratios):
                ax.text(E_over_D, L_over_D, f'{Error_grid[i, j]:.1f}', ha='center', va='center', fontsize=8)

        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label('Dmax Error (%)')

        ax.set_xlabel('E/D Ratio')
        ax.set_ylabel('L/D Ratio')
        ax.set_title(f'Dmax Error Heatmap: {benchmark_analysis} vs {test_analysis}')
        plt.tight_layout()

        if save:
            fig.savefig(save_folder / f'Dmax_Error_Heatmap_{benchmark_analysis}_vs_{test_analysis}.png', dpi=300, bbox_inches='tight')
        if show:
            plt.show()
        plt.close(fig)



######################################################################
        fig, ax = plt.subplots(figsize=(8, 6))

        min_error = np.nanmin(Error_grid)
        max_error = np.nanmax(Error_grid)

        # Safe and unsafe regions
        if min_error < -5 < max_error:
            filled = ax.contourf(E_grid, L_grid, Error_grid, levels=[min_error - 0.001, -5, max_error + 0.001], colors=['lightcoral', 'lightgreen'], alpha=0.8)
            boundary = ax.contour(E_grid, L_grid, Error_grid, levels=[-5], colors='black', linewidths=2)
            ax.clabel(boundary, inline=True, fontsize=9, fmt='%.1f%%')

        elif max_error <= -5:
            ax.set_facecolor('lightcoral')

        else:
            ax.set_facecolor('lightgreen')

        from matplotlib.patches import Patch
        legend_handles = [Patch(facecolor='lightgreen', edgecolor='black', label='Safe: Error ≥ -5%'), Patch(facecolor='lightcoral', edgecolor='black', label='Unsafe: Error < -5%')]
        ax.legend(handles=legend_handles, loc='best')

        ax.set_xlim(np.min(E_over_D_ratios), np.max(E_over_D_ratios))
        ax.set_ylim(np.min(L_over_D_ratios), np.max(L_over_D_ratios))

        ax.set_xlabel('E/D Ratio')
        ax.set_ylabel('L/D Ratio')
        ax.set_title(f'Safe and Unsafe Regions: {benchmark_analysis} vs {test_analysis}')

        plt.tight_layout()

        if save:
            save_folder = Path("Column_Results") / "Analysis_History" / str(frame_id) / "Contour plots"
            save_folder.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_folder / f'Dmax_Error_Safe_Unsafe_{benchmark_analysis}_vs_{test_analysis}.png', dpi=300, bbox_inches='tight')

        if show:
            plt.show()

        plt.close(fig)

#####################################################################

        fig, ax = plt.subplots(figsize=(8, 6))

        min_error = np.nanmin(Error_grid)
        max_error = np.nanmax(Error_grid)

        levels = [min_error]
        if min_error < -5 < max_error or np.isclose(-5, min_error) or np.isclose(-5, max_error):
            levels.append(-5)
        if min_error < 0 < max_error or np.isclose(0, min_error) or np.isclose(0, max_error):
            levels.append(0)
        levels.append(max_error)

        levels = sorted(set(levels))

        if len(levels) >= 2:
            if len(levels) == 2:
                if max_error <= -5:
                    colors = ['lightcoral']
                elif min_error >= 0:
                    colors = ['lightgreen']
                else:
                    colors = ['khaki']
            elif len(levels) == 3:
                if max_error <= 0 and min_error < -5:
                    colors = ['lightcoral', 'khaki']
                elif min_error >= -5 and max_error > 0:
                    colors = ['khaki', 'lightgreen']
                else:
                    colors = ['lightcoral', 'lightgreen']
            else:
                colors = ['lightcoral', 'khaki', 'lightgreen']

            filled = ax.contourf(E_grid, L_grid, Error_grid, levels=levels, colors=colors, alpha=0.8)

        if min_error <= -5 <= max_error:
            contour_minus5 = ax.contour(E_grid, L_grid, Error_grid, levels=[-5], colors='black', linewidths=2)
            ax.clabel(contour_minus5, inline=True, fontsize=9, fmt='%.1f%%')

        if min_error <= 0 <= max_error:
            contour_zero = ax.contour(E_grid, L_grid, Error_grid, levels=[0], colors='blue', linewidths=2)
            ax.clabel(contour_zero, inline=True, fontsize=9, fmt='%.1f%%')

        from matplotlib.patches import Patch
        legend_handles = [
            Patch(facecolor='lightcoral', edgecolor='black', label='Error < -5%'),
            Patch(facecolor='khaki', edgecolor='black', label='-5% ≤ Error < 0%'),
            Patch(facecolor='lightgreen', edgecolor='black', label='Error ≥ 0%')
        ]
        ax.legend(handles=legend_handles, loc='best')

        ax.set_xlabel('E/D Ratio')
        ax.set_ylabel('L/D Ratio')
        ax.set_title(f'Dmax Error Regions: {benchmark_analysis} vs {test_analysis}')
        plt.tight_layout()

        if save:
            save_folder = Path("Column_Results") / "Analysis_History" / str(frame_id) / "Contour plots"
            save_folder.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_folder / f'Dmax_Error_ThreeRegions_{benchmark_analysis}_vs_{test_analysis}.png', dpi=300, bbox_inches='tight')

        if show:
            plt.show()

        plt.close(fig)

def plot_Dmax_error_vs_Radial_Error(json_file_path,benchmark_analysis,test_analysis,L_over_D_ratios, E_over_D_ratios,
                                     frame_names_list=None, storey_heights_list=None, no_of_stories_list=None, bending_axes_list=None, support_list=None,save=True,show=False):

    load_combination_angle = lambda governing_load_combo, L_over_D_ratio, E_over_D_ratio: {
        "1.4D": 90.0,
        "1.2D+1.6L": 90.0,
        "1.2D+L+E": np.degrees(np.arctan2(1.2 + L_over_D_ratio, E_over_D_ratio)),
        "0.9D+E": np.degrees(np.arctan2(0.9, E_over_D_ratio))
    }[governing_load_combo]

    json_files = [f for f in Path(json_file_path).iterdir() if f.is_file() and f.suffix.lower() == ".json"]

    for file_path in tqdm(json_files, desc="Processing Frames", unit="frame"):

        # Find and load json files
        if not file_path.is_file() or file_path.suffix.lower() != ".json":
            continue

        with open(file_path, "r", encoding="utf-8") as file:
            content = json.load(file)

        benchmark_ALR_H = content[benchmark_analysis]["Original_ALR_H"]
        benchmark_ALR_V = content[benchmark_analysis]["Original_ALR_V"]
        benchmark_interaction_diagram = InteractionDiagram2d(benchmark_ALR_H, benchmark_ALR_V)
        test_ALR_H = content[test_analysis]["Original_ALR_H"]
        test_ALR_V = content[test_analysis]["Original_ALR_V"]
        test_interaction_diagram = InteractionDiagram2d(test_ALR_H, test_ALR_V)

        # Filters to run only selected frames
        frame_id = content['data']['Frame_id']
        frame_name = content['data']['column_section']['common_and_exceptions']['common'][0]  ### !!!! This works only for MFColumn Results
        bending_axis = content['data']['column_section']['common_and_exceptions']['common'][1]### !!!! This works only for MFColumn Results
        storey_height = content['data']["storey_height"][0]
        no_of_stories = len(content['data']["storey_height"])
        support = content['data']["support"]

        matches_filter = ((frame_names_list is None or frame_name in frame_names_list) 
                    and (storey_heights_list is None or storey_height in storey_heights_list)
                    and (no_of_stories_list is None or no_of_stories in no_of_stories_list) 
                    and (bending_axes_list is None or bending_axis in bending_axes_list) 
                    and (support_list is None or support in support_list))

        if not matches_filter:
            continue


        for L_over_D in tqdm(L_over_D_ratios, desc=f"Frame {frame_id}", unit="ratio", leave=False):

            results_benchmark = get_Dmax_EoverD_LoverD(content, benchmark_analysis, L_over_D_ratios=[L_over_D], E_over_D_ratios=E_over_D_ratios)
            results_test = get_Dmax_EoverD_LoverD(content, test_analysis, L_over_D_ratios=[L_over_D], E_over_D_ratios=E_over_D_ratios)
            E_over_D = [v['E_over_D_ratio'] for v in results_benchmark.values()]
            Dmax_benchmark = [v['Dmax'] for v in results_benchmark.values()]
            Dmax_test = [v['Dmax'] for v in results_test.values()]
            Governing_combo_test = [v['Governing_Load_Combo'] for v in results_test.values()]
            angles = [load_combination_angle(combo, L_over_D, E) for combo, E in zip(Governing_combo_test, E_over_D)]
            radial_error = benchmark_interaction_diagram.compare_two(test_interaction_diagram, angles=angles, degrees=True)
            Dmax_error = [(benchmark - test) / benchmark for benchmark, test in zip(Dmax_benchmark, Dmax_test)]
            plt.scatter(radial_error, Dmax_error, s=15)

        plt.xlabel("Radial Error (%)")
        plt.ylabel("Dmax Error (%)")
        plt.grid(True)
        plt.show()

def compute_del2_over_del1(json_file_path,theta_list,analysis_types_list=None):
    from MFColumn2D import MFColumn_2D
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


        analyses = [
            key
            for key, value in content.items()
            if key != "data" and isinstance(value, dict)
        ]

        if not analyses:
            print(f"Skipping {file_path.name}: no analyses were found.")
            continue

        for analysis_index, analysis in enumerate(analyses):
            matches_filter = (analysis_types_list is None or analysis in analysis_types_list)

            if not matches_filter:
                continue
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
                Material_dict=Material_Info[data['Steel_Grade']]
                Material_details=convert_dict_items_to_class_attributes(Material_dict)
                Steel=Steel_Material(mat_tag=1,E=Material_details.E,Fy=Material_details.Fy)
                data['Material_obj']= Steel
                Frame = MFColumn_2D.from_dict(data)
                print(Frame.P_M_M_interaction_limit)
                del2_over_del1=Frame.get_del2_over_del1(vertical_load_scale=vertical_load_scale, lateral_load_scale=lateral_load_scale)
                results['Theta'].append(theta)
                results['ALR_H'].append(pt[0] if pt is not None else None)
                results['ALR_V'].append(pt[1] if pt is not None else None)
                results['del2_over_del1'].append(del2_over_del1)

            updated_analysis_data = analysis_data.copy()
            updated_analysis_data["del2_over_del1_results"] = results

            Results = {analysis: updated_analysis_data}
            update_or_create_json(file_path, Results)

            
def plot_error_vs_del2_over_del1(json_file_path, analysis_types, cap=3,
                                 frame_names_list=None, storey_heights_list=None,
                                 no_of_stories_list=None, bending_axes_list=None,
                                 support_list=None, save=True, show=True):

    Folder = Path(json_file_path)
    save_folder = os.path.join("Column_Results", "error_vs_del2_over_del1")
    os.makedirs(save_folder, exist_ok=True)

    benchmark_analysis = analysis_types[0]

    for test_analysis in analysis_types[1:]:

        plt.figure(figsize=(8, 6))

        for file_path in tqdm(sorted(Folder.glob("*.json")), desc=f"Processing {test_analysis}", unit="frame"):

            with open(file_path, "r", encoding="utf-8") as file:
                content = json.load(file)

            # Frame details
            data = content['data']
            frame_name = data['column_section']['common_and_exceptions']['common'][0]
            bending_axis = data['column_section']['common_and_exceptions']['common'][1]
            storey_height = data['storey_height'][0]
            no_of_stories = len(data['storey_height'])
            support = data['support']

            # Apply filters
            matches_filter = (
                (frame_names_list is None or frame_name in frame_names_list)
                and (storey_heights_list is None or storey_height in storey_heights_list)
                and (no_of_stories_list is None or no_of_stories in no_of_stories_list)
                and (bending_axes_list is None or bending_axis in bending_axes_list)
                and (support_list is None or support in support_list)
            )

            if not matches_filter:
                continue

            # Slenderness ratio
            r = (wide_flange_database[frame_name]['rx'] if bending_axis == 'x' else wide_flange_database[frame_name]['ry']) * inch
            slenderness_ratio = storey_height / r

            # Interaction diagrams
            benchmark = InteractionDiagram2d(content[benchmark_analysis]['Original_ALR_H'], content[benchmark_analysis]['Original_ALR_V'])
            test = InteractionDiagram2d(content[test_analysis]['Original_ALR_H'], content[test_analysis]['Original_ALR_V'])

            # Calculate errors
            theta = content[test_analysis]['del2_over_del1_results']['Theta']
            del2_over_del1 = np.array(content[test_analysis]['del2_over_del1_results']['del2_over_del1'], dtype=float)
            errors = np.array(benchmark.compare_two(test, theta, degrees=True), dtype=float)

            # Filter invalid data and values exceeding cap
            valid = np.isfinite(del2_over_del1) & np.isfinite(errors) & (del2_over_del1 <= cap)
            X, Y = del2_over_del1[valid], errors[valid]

            # Sort data for continuous curves
            order = np.argsort(X)
            color = plt.cm.viridis((slenderness_ratio - 10) / (120 - 10))
            plt.plot(X[order], Y[order], color=color, linewidth=1.8)

        # Plot formatting
        plt.axhline(0, color='black', linewidth=1)
        plt.axhline(-0.05, color='red', linestyle='--', linewidth=1)
        plt.xlabel(r'$\Delta_2/\Delta_1$')
        plt.ylabel('Radial Error')
        plt.title(f'{benchmark_analysis} vs {test_analysis}')
        sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=10, vmax=120))
        plt.colorbar(sm, ax=plt.gca(), label=r'Slenderness Ratio, $L_c/r$')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        if save:
            plt.savefig(os.path.join(save_folder, f"Error_vs_del2_over_del1_{benchmark_analysis}_vs_{test_analysis}.png"), dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        else:
            plt.close()


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


def save_analysis_history_figures(parent_folder, results, frame_name, analysis_type, ALR_H, ALR_V,
                                  step_size=None, tolerance=None, iterations=None, control_dir=None):

    parent_folder = Path(parent_folder)
    parent_folder.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(16, 9))
    plot_options = {"linewidth": 2, "marker": "o", "markersize": 3, "markevery": 1}
    # Plot definitions: (x, y, xlabel, ylabel, title)
    plots = [
        (results.control_node_displacement, results.load_ratio,
         "Displacement at Control Node", "Load Ratio λ", "Load Ratio vs Displacement"),
        (results.load_ratio, results.lowest_eigenvalue,
         "Load Ratio λ", "Lowest Eigenvalue", "Eigenvalue vs Load Ratio"),
        (results.vertical_reaction, results.load_ratio,
         "Vertical Reaction", "Load Ratio λ", "Load Ratio vs Vertical Reaction"),
        (results.base_shear, results.load_ratio,
         "Base Shear", "Load Ratio λ", "Load Ratio vs Base Shear"),
        (results.absolute_maximum_strain, results.load_ratio,
         "Absolute Maximum Strain", "Load Ratio λ", "Load Ratio vs Maximum Strain"),
        (results.load_ratio, results.max_P_M_M_interaction,
         "Load Ratio λ", "Maximum P-M-M Interaction", "P-M-M Interaction vs Load Ratio")
    ]
    # Plot all six graphs
    for i, (x, y, xlabel, ylabel, title) in enumerate(plots, start=1):
        plt.subplot(2, 3, i)
        plt.plot(x, y, **plot_options)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True, linestyle="--", linewidth=0.7, alpha=0.5)
        plt.tick_params(axis="both", direction="in", labelsize=9)
        plt.margins(x=0.05, y=0.08)
    # Main title
    title = f"{frame_name} — {analysis_type}"
    if control_dir is not None:
        title += f" — Control Direction: {control_dir}"

    plt.suptitle(title, fontsize=15, fontweight="bold", y=0.98)
    # Information box
    information = [f"ALR_H: {ALR_H}", f"ALR_V: {ALR_V}"]
    if tolerance is not None:
        information.append(f"Tolerance: {tolerance:.2e}")
    if iterations is not None:
        information.append(f"Iterations: {iterations}")
    if step_size is not None:
        information.append(f"No_of_steps: {step_size}")
    information_text = "\n".join(information)
    plt.gcf().text(0.985, 0.975, information_text, ha="right", va="top", fontsize=9,
                   bbox={"boxstyle": "round,pad=0.4", "facecolor": "white",
                         "edgecolor": "black", "linewidth": 0.8, "alpha": 0.9})
    plt.tight_layout(rect=[0, 0, 1, 0.91], h_pad=2, w_pad=2)
    # Dynamic filename
    filename_parts = [f"ALR_H_{ALR_H}", f"ALR_V_{ALR_V}"]
    if control_dir is not None:
        filename_parts.append(f"control_{control_dir}")
    if tolerance is not None:
        filename_parts.append(f"tol_{tolerance:.2e}".replace("+", ""))
    if iterations is not None:
        filename_parts.append(f"iter_{iterations}")
    if step_size is not None:
        filename_parts.append(f"No_of_steps_{step_size}")
    filename = parent_folder / ("_".join(filename_parts) + ".png")
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()

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



if __name__ == "__main__":
    print("running")
    json_file_path="Column_Results/json_files"
    L_over_D_ratios=np.arange(0,6,0.1)
    # L_over_D_ratios=[1]
    E_over_D_ratios=np.arange(0,2,0.1)
    # E_over_D_ratios=[2]
    analysis_types=['GMNIA','GNA','GNA_Notional_Loads']

    # for analysis_type in analysis_types:
    #     get_governing_load_combo_details(json_file_path,analysis_type,L_over_D_ratios,E_over_D_ratios)
    
    # plot_interaction_diagrams_for_all_analyses_with_load_combo_info(json_file_path=json_file_path,
    #                                                                 use_original=True,
    #                                                                 save=True,
    #                                                                 show=False,
    #                                                                 analysis_types=analysis_types,
    #                                                                 load_combo_analysis='GNA',
    #                                                                 load_combo_L_over_D_ratio=None,
    #                                                                 load_combo_E_over_D_ratio=None)  


    #plot_governing_D_max_vs_E_over_D_for_all_analyses(json_file_path=json_file_path,save=True,show=False,analysis_types=analysis_types,L_over_D_ratios=L_over_D_ratios,E_over_D_ratios=E_over_D_ratios)

    plot_LoverD_EoverD_Dmaxerror_contour_plots(json_file_path,'GMNIA','GNA',L_over_D_ratios,E_over_D_ratios)
    # plot_Dmax_error_vs_Radial_Error(json_file_path,'GMNIA','GNA',L_over_D_ratios,E_over_D_ratios)