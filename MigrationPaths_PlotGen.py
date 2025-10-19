# Sami Alawadhi, David M. Rutkowski (2025) 
# Vavylonis Group
# Department of Physics, Lehigh University
# If you use any part of the scripts in this package, please cite:
# TBD

import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
import pandas as pd  # Import pandas for writing to Excel

def read_migration_data(filepath, num_NuclearBeads):
    with open(filepath) as fp:
        XYZ_line = fp.readline()
        count = 0
        x = []
        y = []
        while XYZ_line:
            split_XYZ_string = XYZ_line.split() 

            if len(split_XYZ_string) == 1:
                num_beads = int(split_XYZ_string[0])
                XYZ_line = fp.readline()

                sum_cmY = 0.0
                sum_cmX = 0.0 
                for b in range(num_beads):
                    XYZ_line = fp.readline()
                    split_XYZ_string = XYZ_line.split()
                    if int(split_XYZ_string[0]) == 0:
                        sum_cmY += float(split_XYZ_string[2])
                        sum_cmX += float(split_XYZ_string[3])
                x.append(-sum_cmX / num_NuclearBeads)
                y.append(sum_cmY / num_NuclearBeads)
                count += 1
            XYZ_line = fp.readline()
        x = np.array(x) - x[0]  # Normalize x to start at 0
        y = np.array(y) - y[0]  # Normalize y to start at 0
        return x, y

def setup_file_selection(multiple=True):
    root = tk.Tk()
    root.withdraw()
    filepaths = []
    while True:
        filepath = filedialog.askopenfilename()
        if filepath == "":
            break
        filepaths.append(filepath)
    root.destroy()
    return filepaths

def save_to_excel(filepaths, writer, compartmentalization):
    num_NuclearBeads = 100.0  # Assuming 100 beads for nucleus
    for i, filepath in enumerate(filepaths, start=1):
        x, y = read_migration_data(filepath, num_NuclearBeads)
        # Create a DataFrame with x and y values
        df = pd.DataFrame({'X (µm)': x, 'Y (µm)': y})
        sheet_name = f"Grid{i}"
        df.to_excel(writer, sheet_name=sheet_name, index=False)
        print(f"Saved data for {compartmentalization} Grid {i} to sheet '{sheet_name}'")

def plot_paths(filepaths, color, label):
    num_NuclearBeads = 100.0  
    for filepath in filepaths:
        x, y = read_migration_data(filepath, num_NuclearBeads)
        plt.plot(x, y, color=color, label=label if filepath == filepaths[0] else "")

def plot_migration_paths():
    # File selection
    print("Select the files for septin compartmentalization, press Cancel when done:")
    septin_filepaths = setup_file_selection(multiple=True)
    print("Select the files for no compartmentalization, press Cancel when done:")
    no_septin_filepaths = setup_file_selection(multiple=True)

    # Color-blind friendly colors
    colors = {
        "CB": "blue",  
        "NoCB": "red"
    }

    # Save data to separate Excel files
    with pd.ExcelWriter("Migration_Paths_CB.xlsx", engine="xlsxwriter") as writer_cb:
        save_to_excel(septin_filepaths, writer_cb, "CB")

    with pd.ExcelWriter("Migration_Paths_NoCB.xlsx", engine="xlsxwriter") as writer_nocb:
        save_to_excel(no_septin_filepaths, writer_nocb, "NoCB")

    # Data processing and plotting:
    plot_paths(septin_filepaths, colors["CB"], '- Compartmentalized T-cells')
    plot_paths(no_septin_filepaths, colors["NoCB"], '- Non-compartmentalized T-cells')

    #plt.scatter([0], [0], color='black', s=30, zorder=5, label='Start Position')  # Start position marker
    plt.plot(0, 0, marker='+', color='darkred', markersize=12, markeredgewidth=5, zorder=5)  # Start position marker

    # Orange dashed cross lines stopping before the origin
    plt.axhline(0, color='orange', linestyle='--', linewidth=2, xmin=0.1, xmax=0.9)
    plt.axvline(0, color='orange', linestyle='--', linewidth=2, ymin=0.05, ymax=0.95)

    plt.title(r'$\bf{+/- compartmentalization}$, No repolarization', fontsize=16, pad=10)
    plt.xlabel('Distance (\u03BCm)', fontsize=14)
    plt.tick_params(axis='x', labelsize=12, direction='in', pad=-20)  # Move x-axis numbers inside
    plt.ylabel('Distance (\u03BCm)', fontsize=14)
    plt.tick_params(axis='y', labelsize=12, direction='in', pad=-30)  # Move y-axis numbers inside
    plt.gca().text(155, 0, 'Distance (\u03BCm)', fontsize=14, rotation=-90, va='center')
    plt.legend(loc='lower right', fontsize=12, frameon=True, facecolor='white', edgecolor='lightgrey')
    plt.grid(True, linestyle='-', linewidth=0.3)

    plt.xlim(-150, 150)
    plt.ylim(-150, 150)
    plt.xticks(np.arange(-100, 101, 50))  # Set x-axis ticks from -100 to 100 with a step of 50
    plt.yticks(np.arange(-100, 101, 50))  # Set y-axis ticks from -100 to 100 with a step of 50

    plt.gca().set_aspect('equal', adjustable='box')
    plt.gca().tick_params(axis='both', direction='in', top=True, right=True, labeltop=True, labelright=True)
    plt.figure(figsize=(6, 6), dpi=600)

    plt.show()

# Main function to run
plot_migration_paths()
