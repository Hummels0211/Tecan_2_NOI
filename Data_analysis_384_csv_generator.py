import pandas as pd
from pathlib import Path
import re
from openpyxl import Workbook, worksheet, load_workbook
import os

# Set the current working directory
os.chdir("/home/jiayu/Documents/Metabolites_info/Single_met_pref_test/Tecan_384/Exp_results_B2")

# Load the mapping files
mapping_norm = pd.read_csv("./Rand_mapping/Randomised_384_full_plate_norm.csv", sep = ",")
mapping_mod = pd.read_csv("./Rand_mapping/Randomised_384_full_plate_mod.csv", sep = ",")

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

# Load the dynamics from each well by looping through the .xlsx file
#data_b1 = pd.read_excel('./OD_measurement/OD_600_C1_C4_10_02_2023.xlsx', sheet_name='Sheet2', skiprows=58, nrows=300)
#data_b2 = pd.read_excel('./OD_measurement/OD_600_C5_C8_14_03_2023.xlsx', sheet_name='Sheet2', skiprows=58, nrows=300)
#data_b3 = pd.read_excel('./OD_measurement/OD_600_C9_C12_10_08_2023.xlsx', sheet_name='Sheet2', skiprows=58, nrows=300)
#data_b4 = pd.read_excel('./OD_measurement/OD_600_C13_C16_17_08_2023.xlsx', sheet_name='Sheet2', skiprows=58, nrows=300)
#data_b5 = pd.read_excel('./OD_measurement/OD_600_C17_C20_19_08_2023.xlsx', sheet_name='Sheet2', skiprows=58, nrows=300)

# Create a dictionary to store all the data as dataframes

df_data_b1 = pd.DataFrame()
df_data_b2 = pd.DataFrame()
df_data_b3 = pd.DataFrame()
df_data_b4 = pd.DataFrame()
df_data_b5 = pd.DataFrame()

df_dict = [df_data_b1, df_data_b2, df_data_b3, df_data_b4, df_data_b5]

# Generate a final reading table by looping the excel file
for (root, dirs, files) in os.walk("./OD_measurement/"):
    for xlsx, j in zip(files, range(5)):
        filepath = os.path.join(root, xlsx)
        for i in range(384):
            df = pd.read_excel(filepath, sheet_name='Sheet2', skiprows=(58 + 303 * i), nrows=300)
            header = list(df)
            # Change the first column as the well labels
            df.iloc[:, 0] = header[0]
            df = df.rename(columns={header[0]: 'Well_label'})
            # Remove the unecessary columns
            df = df.drop(columns=['Temp. [°C]', '0;1', '1;1', '1;0', '0;0'])
            # Expand the final dataframe
            df_dict[j] = pd.concat([df_dict[j], df], ignore_index=True)
        # Save the dataframe to a new csv file
        df_dict[j].to_csv(f"./batch_{j+1}.csv", sep = ",", index = False)
