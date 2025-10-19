# Sami Alawadhi, David M. Rutkowski (2025) 
# Vavylonis Group
# Department of Physics, Lehigh University
# If you use any part of the scripts in this package, please cite:
# TBD

import numpy as np
from matplotlib import pyplot as plt
import tkinter as tk
from tkinter import filedialog
import pandas as pd  # Import pandas for Excel writing

# Set to keep track of labels already used in the legend
used_labels = set()

def plot_trajectory(XYZ_filepath, color, linestyle, compartmentalization, grid_num, writer):
    global used_labels  # Use the global set to manage labels

    # Initialize the obstacle grid entry point
    obstacle_entry_point = -17.39
    entry_index = None

    with open(XYZ_filepath) as fp:
        XYZ_line = fp.readline()
        count = 0
        snapshotTime = 30.84 # seconds
        scaleTime = 1.0
        t = []
        x = []
        num_NuclearBeads = 100.0
        while XYZ_line:
            split_XYZ_string = XYZ_line.split()  # Split elements on different spacing.

            if len(split_XYZ_string) == 1:
                num_beads = int(split_XYZ_string[0])
                XYZ_line = fp.readline()

                sum_x = 0.0
                for b in range(num_beads):
                    XYZ_line = fp.readline()
                    split_XYZ_string = XYZ_line.split()
                    if int(split_XYZ_string[0]) == 0:
                        sum_x += float(split_XYZ_string[3])

                cm_x = -(sum_x) / num_NuclearBeads
                if entry_index is None and cm_x < obstacle_entry_point:
                    entry_index = count
                    entry_cm_x = cm_x

                t.append(((count * snapshotTime)/scaleTime) / 60)  # Time in minutes
                x.append(cm_x)
                count += 1

                if count >= 8000:  # Stop counting once 8000 is reached- modify as needed.
                    break

            XYZ_line = fp.readline()

    x = [i - x[0] for i in x]  # Adjust displacement to start from zero

    if entry_index is not None:
        t = [time - t[entry_index] for time in t]
        x = [position - entry_cm_x for position in x]

    # Print the (time, displacement) pairs
    time_displacement_pairs = list(zip(t, x))
    print(f"Time and Displacement values for {compartmentalization} Grid {grid_num}: {time_displacement_pairs}")

    # Save the time and displacement values to a subsheet in the Excel file
    if compartmentalization == "With compartmentalization":
        sheet_name = f"Grid{grid_num}_CB"  # CB for Compartment Boundaries
    else:
        sheet_name = f"Grid{grid_num}_NoCB"  # NoCB for No Compartment Boundaries

    df = pd.DataFrame(time_displacement_pairs, columns=["Time (min)", "Displacement (µm)"])
    df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"Saved {compartmentalization} data for Grid {grid_num} to subsheet '{sheet_name}'")

    # Label handling: only add if not already used
    label = compartmentalization if compartmentalization not in used_labels else ""

    plt.plot(t, x, color=color, linestyle=linestyle, label=label)
    used_labels.add(compartmentalization)  # Mark this label as used

# Color-blind friendly colors
colors = {
    "With compartmentalization": "#0006B2",  # Blue
    "Without compartmentalization": "#D50000"  # Red
}

root = tk.Tk()
root.withdraw()

# Create Excel writers for both conditions
with pd.ExcelWriter("With_Compartmentalization_rescaledTime.xlsx", engine="xlsxwriter") as with_writer, \
     pd.ExcelWriter("Without_Compartmentalization_rescaledTime.xlsx", engine="xlsxwriter") as without_writer:

    for grid_num in range(1, 13):  # Loop over N-1 grids
        # For Septin Compartmentalization
        XYZ_filepath = filedialog.askopenfilename(title=f"Select file for Grid {grid_num}, With compartmentalization")
        if XYZ_filepath != "":
            plot_trajectory(XYZ_filepath, color=colors["With compartmentalization"], linestyle="-", compartmentalization="With compartmentalization", grid_num=grid_num, writer=with_writer)

        # For No Compartmentalization
        XYZ_filepath = filedialog.askopenfilename(title=f"Select file for Grid {grid_num}, Without compartmentalization")
        if XYZ_filepath != "":
            plot_trajectory(XYZ_filepath, color=colors["Without compartmentalization"], linestyle="-", compartmentalization="Without compartmentalization", grid_num=grid_num, writer=without_writer)

# Set the axis limits to ensure (0,0) is at the corner
plt.xlim(left=0)
plt.ylim(bottom=0)

plt.gca().tick_params(axis='both', labelsize=18)

plt.xlabel('Time (min)', fontsize=20)
plt.tick_params(axis='x', labelsize=20)
plt.ylabel('Displacement (\u03BCm)', fontsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.legend(loc='best', fontsize=18)
plt.show()
