import pandas as pd
from pathlib import Path
import re
import os

# To generate the plots
import matplotlib.pyplot as plt
import seaborn as sns

# Set the default working directory
os.chdir("/home/jiayu/Documents/Metabolites_info/Single_met_pref_test/Tecan_384/Exp_results_B2")

# Define a function to match the well numbers (1-384) to the well labels (A1-P24)
def wellnumber2label(x):
    col_num = (x - 1)//16 + 1
    row_lab = chr((x - 1)%16 + 65) # Capital Letter 'A' starts from 65 in ASCII system
    well = row_lab + str(col_num)
    return well

# Define a function convert well labels (A1-P24) to well numbers (1-384)
def welllabel2number(x):
    split = re.split('(\d+)', x)
    row = ord(split[0]) - 64
    col = int(split[1]) - 1
    well = col * 16 + row
    return well

# Load the mapping files
#mapping_norm = pd.read_csv("./Rand_mapping/Randomised_384_full_plate_norm.csv", sep = ",")
#mapping_mod = pd.read_csv("./Rand_mapping/Randomised_384_full_plate_mod.csv", sep = ",")

# Load the measurements data
df_data_b1 = pd.read_csv('./batch_1.csv', sep = ',')
df_data_b2 = pd.read_csv('./batch_2.csv', sep = ',')
df_data_b3 = pd.read_csv('./batch_3.csv', sep = ',')
df_data_b4 = pd.read_csv('./batch_4.csv', sep = ',')
df_data_b5 = pd.read_csv('./batch_5.csv', sep = ',')

# Assign the mapping files 
mapping_b1 = pd.read_csv("./Rand_mapping/Randomised_384_full_plate_norm.csv", sep = ",")
mapping_b2 = pd.read_csv("./Rand_mapping/Randomised_384_full_plate_norm.csv", sep = ",")
mapping_b3 = pd.read_csv("./Rand_mapping/Randomised_384_full_plate_norm.csv", sep = ",")
mapping_b4 = pd.read_csv("./Rand_mapping/Randomised_384_full_plate_norm.csv", sep = ",")
mapping_b5 = pd.read_csv("./Rand_mapping/Randomised_384_full_plate_mod.csv", sep = ",")

# Test the programme
#print(mapping_b2)

mapping_list = [mapping_b1, mapping_b2, mapping_b3, mapping_b4, mapping_b5]
data_list = [df_data_b1, df_data_b2, df_data_b3, df_data_b4, df_data_b5]

# Modify the numbers of the metabolites for carbon list mapping
for i in range(len(mapping_list)):
    mapping_list[i]['Carbon_number'] = mapping_list[i]['Carbon_number'] + (i * 4)
    # Modify the label of the wells for merging with the measurment data
    mapping_list[i]['Well_label'] = mapping_list[i]['Well_number'].apply(lambda x: wellnumber2label(x))

# Test the programme
#print(len(mapping_list))

# Merging the data and the mapping list
for i in range(len(mapping_list)):
    # Access each metadatafile from data_list variable
    data_list[i] = pd.merge(data_list[i], mapping_list[i], how = 'left', on = 'Well_label')


# Put all the measurements together as a metadata table
df_meta = pd.concat(data_list)

# Load the mapping file for the SynCom strains
sc_mapping = pd.read_csv('/home/jiayu/Documents/SynCom_Info/Syncom_mapping_2019.csv', sep = ',', index_col = 0)
sc_mapping = sc_mapping[['ID', 'phylum', 'class', 'order', 'family', 'genus']]
sc_mapping = sc_mapping.rename({'ID': 'Bacteria_Strain'}, axis=1)

# Load the treatment mapping list
metabolite_mapping = pd.read_csv('/home/jiayu/Documents/Metabolites_info/Single_met_pref_test/Tecan_384/Exp_design/metabolite_mapping.tsv', sep = '\t')

# Test the programme
#print(mapping_b2)

# Merge the data and taxanomy and treatment mapping lists
df_meta = pd.merge(df_meta, sc_mapping, how = 'left', on = 'Bacteria_Strain')
df_meta = pd.merge(df_meta, metabolite_mapping, how = 'left', on = 'Carbon_number')
df_meta['Time [h]'] = df_meta['Time [s]'] / 3600

# Save the merged table to a new metadata csv file
df_meta.to_csv('./meta_table.csv', sep = ',', index=False)

# Generate all the growth curves for each combination
if not os.path.isdir('./Growth_curve_plots'):
   os.makedirs('./Growth_curve_plots')

# Setting conditions for looping the metadata
for i in df_meta['Bacteria_Strain'].unique():
        for j in df_meta['Treatment'].unique():
            cond_1 = df_meta['Bacteria_Strain'] == i
            cond_2 = df_meta['Treatment'] == j
            df = df_meta[cond_1 & cond_2]
            fig = plt.figure(figsize=(8,6))
            ax = plt.axes()
            sns_plot = sns.lineplot(data = df, x = 'Time [h]', y = 'Mean', hue='Well_label')
            plt.title(i + '_' + j)
            plt.tight_layout()
            plt.savefig('./Growth_curve_plots/' + i + '_' + j + '_dyn.png', dpi=300)
