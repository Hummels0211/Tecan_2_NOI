import pandas as pd
import seaborn as sns
import scipy as sp
import numpy as np
import os
import re
import matplotlib.pyplot as plt
from pathlib import Path

# Set the default working directory
os.chdir("/home/jiayu/Documents/Metabolites_info/Single_met_pref_test/Tecan_384/Exp_results_B2/")

# Load the metadata table containing all the results
df_result = pd.read_csv("./meta_table.csv", sep = ",", header = 0)

# Load the SynCom mapping file
df_mapping = pd.read_csv("~/Documents/SynCom_Info/Syncom_mapping_2019.csv", sep = ",")

# Create an empty dataframe to store the first derivative information
df_result_gr = pd.DataFrame()

# Calibrate the OD measurements by substrating the OD value at T_0
for ind, row in df_result.iterrows():
    i = (ind // 300) * 300
    df_result.loc[ind, "Norm_Mean"] = df_result.loc[ind, "Mean"] - df_result.loc[i, "Mean"]

# Calculate the real-time growth rate 
for well in df_result["Well_label"].unique():
    df_temp = df_result[df_result["Well_label"] == well]
    for carbon in df_temp["Treatment"].unique():
        df_temp_2 = df_result.loc[(df_result["Well_label"] == well) & (df_result["Treatment"] == carbon)]
        df_temp_2["Growth_rate"] = df_temp_2["Mean"].diff()
        df_result_gr = pd.concat([df_result_gr, df_temp_2])

# Save the results with growth rates to a new .csv file to test
df_result_gr.to_csv("./meta_table_gr.csv", sep = ",")

# pair = df_result[["Well_label", "Treatment"]].drop_duplicates()

# Generate all the growth curves for each combination
if not os.path.isdir('./Snapshot_data'):
   os.makedirs('./Snapshot_data')

# Set the desired numbers of snapshots from the full data
snapshots = 5

# Generate .csv files using loops to capture the growth every n hours
for i in range(-1, 300, 300//snapshots):
    hr = i // 6 + 1
    if i < 0:
        i = 0
        df = df_result_gr.iloc[i:len(df_result_gr):300, :]
        # Use the derivative at the first time point to avoid the boundary scenes
        df_2 = df_result_gr.iloc[i+1:len(df_result_gr):300, :]
    else:
        df = df_result_gr.iloc[i:len(df_result_gr):300, :]
        df_2 = df_result_gr.iloc[i:len(df_result_gr):300, :]
    df = df.drop(df[df.Bacteria_Strain == "Pos_Ctrl"].index)
    df_2 = df_2.drop(df_2[df_2.Bacteria_Strain == "Pos_Ctrl"].index)
    df_gr = df_2[["Well_label", "Bacteria_Strain", "Treatment", "Growth_rate"]]
    df = df[["Well_label", "Bacteria_Strain", "Treatment", "Norm_Mean"]]
    df = df.groupby(["Well_label", "Bacteria_Strain", "Treatment"]).mean()
    df_gr = df_gr.groupby(["Well_label", "Bacteria_Strain", "Treatment"]).mean()
    df = df.reset_index(level = [0,1,2])
    df_gr = df_gr.reset_index(level = [0,1,2])
    # Calculate the integral of the growth curves by using pandas.agg()
    df = df.groupby(["Bacteria_Strain", "Treatment"]).agg({'Norm_Mean':['mean','std']})
    df.columns = df.columns.get_level_values(1)
    df = df.reset_index()
    df.to_csv('./Snapshot_data/OD_snapshot_{}h.csv'.format(hr), index = False)
    df_gr.to_csv('./Snapshot_data/Growth_rate_snapshot_{}h.csv'.format(hr), index = False)
    

# Nature sorting on the values (strain names) from the columns from the dataframe
import natsort as ns

if not os.path.isdir("./Heatmaps"):
    os.makedirs("./Heatmaps")

import glob

# Create the mapping 
map_family = dict(df_mapping.set_index(["ID", "family"]).index.unique())
# Remove Chlamy165 due to the duplication with Chlamy123
del map_family['Chlamy165']
# Create the mapping 
family_palette = sns.hls_palette(len(df_mapping['family'].unique()), s = 0.35)
family_lut = dict(zip(map(str, df_mapping['family'].unique()), family_palette))
family_colour = pd.Series(map_family).map(family_lut)

# Generate heatmaps by looping over the snapshot datasets
for filename in glob.iglob("./Snapshot_data/OD_snapshot_*h.csv"):
    # sub_dir = "./Snapshot_data/"
    dig_pattern = re.compile(r"\d+")
    time_point = dig_pattern.findall(filename)[0]
    df = pd.read_csv(filename, sep = ",")
    # Resize the dataframe
    df = df.reset_index().pivot(columns = "Bacteria_Strain", index = "Treatment", values = "mean")
    # Define the plot
    fig, ax = plt.subplots(figsize=(13,7))
    # Add title to the Heat map
    title = "Metabolite Preference Map" + time_point + "h"
    # Set the font size and the distance of the title from the plot
    plt.title(title,fontsize=18)
    ttl = ax.title
    ttl.set_position([0.5,1.05])
    # Hide ticks for X & Y axis
    #ax.set_xticks([])
    #ax.set_yticks([])
    # Remove the axes
    #ax.axis('off')
    sns.heatmap(data = df, cmap="BuPu", linewidths=2, linecolor='white', ax = ax, cbar_kws = {"shrink": 0.75},
            vmin = 0.0, vmax = 0.7)
    plt.xticks(rotation=45)
    plt.savefig("./Heatmaps/heatmap_" + time_point + "h.png",  dpi = 400)
    # creating mask
    mask = np.triu(np.ones_like(df.corr(method = 'spearman')))
    # Define the triangle plot of correlation
    fig2, ax2 = plt.subplots(figsize = (9,8))
    sns.heatmap(df.corr(method = "spearman"), linewidth = 2.5, linecolor = 'white', square = True, #
            cbar_kws = {"shrink": 0.6, "label": "Spearman Correlation"}, cmap = "BuGn", mask = mask,
            vmin=0, vmax=1)
    plt.title('Pair-wise correlation of metabolite preference'+" ("+time_point+"h)", fontsize = 18)
    plt.xlabel('Bacteria Strain', fontsize = 15)
    plt.ylabel('Bacteria Strain', fontsize = 15)
    plt.xticks(rotation=45, ha = 'right')
    plt.yticks(rotation=45, va = 'top')
    plt.tight_layout()
    plt.savefig("./Heatmaps/correlation_" + time_point + "h.png", dpi = 400)
    # Plot the clustered heatmap with taxonomy table
    fig_3 = sns.clustermap(df.corr(method='spearman'), linewidth = 2.5, linecolor = 'white', square = True,
                     cbar_kws = {"shrink": 0.1, "orientation": "horizontal", "ticks": [-1, 0, 1]},
                     cmap = "RdYlGn", col_colors= family_colour,
                     vmin=-1, vmax=1)
    x0, _y0, _w, _h = fig_3.cbar_pos
    fig_3.ax_cbar.set_position([x0, 0.95, fig_3.ax_row_dendrogram.get_position().width, 0.02])
    fig_3.ax_cbar.set_title('Spearman Coefficient at ' + time_point + 'h')
    fig_3.ax_cbar.tick_params(axis='x', length=10)
    for spine in fig_3.ax_cbar.spines:
        fig_3.ax_cbar.spines[spine].set_color('#d98880')
        fig_3.ax_cbar.spines[spine].set_linewidth(1)
    for label in df_mapping['family'].unique():
        fig_3.ax_row_dendrogram.bar(0, 0, color=family_lut[label], label=label, linewidth=0)
    fig_3.ax_row_dendrogram.legend(loc="center", ncol=1, bbox_to_anchor=(6.5, 1.), fontsize = 12)
    plt.setp(fig_3.ax_heatmap.get_xticklabels("Bacteria Strain"), rotation=45, ha = 'right')
    plt.savefig("./Heatmaps/correlation_clustered_" + time_point + "h.png", dpi = 500)
