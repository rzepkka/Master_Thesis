import pandas as pd
import numpy as np
import pickle
import json

# ==========================================================================================================================================================================
# =====================================================================================================================================================================
# =========== ADD ALL INPUT DATA FILES HERE ======================================================================================================
# =====================================================================================================================================================================
# =====================================================================================================================================================================

# =================== 1. LOAD PICKLE FILE ========================================================================================================================================================

# add path to your pickle file
read_input_file = open('data/Data.pickle','rb')
load_inputs = pickle.load(read_input_file)
read_input_file.close()

T, S, Sboot = load_inputs

# =================== 2. LOAD PATIENT DATA ========================================================================================================================================================

# add path to your csv file with patients' data
data = pd.read_csv("data/EDADS_data.csv")

# load diagnosis as numpy array 
diagnosis = np.array(data['Diagnosis'])

# =================== 3. ADD DISEASE NAME AND SUBTYPE LABELS =====================================================================================================================================

# pass disease name
disease = 'AD'

# choose subtype names
subtype_labels = ['Subcortical subtype', 'Frontal subtype', 'Parietal subtype','Typical subtype']

# =================== 4. LOAD JSON FILES FOR BRAIN REGION MAPPINGS ===============================================================================================

# 3D mappings
f = open('data/DK_3D_combined.json')
map_dk = json.load(f)
f.close()

f = open('data/ASEG_3D_combined.json')
map_aseg = json.load(f)
f.close()

# 2D mappings
f = open('data/DK_2D_combined.json')
map_dk_2D = json.load(f)
f.close()

f = open('data/ASEG_2D_combined.json')
map_aseg_2D = json.load(f)
f.close()


# ======= 5. LOAD PATH TO PDF FILE =============================================================================================================================================

# add path to your pdf
path_to_pdf = "/Users/macos/Documents/GitHub/Master_Thesis/data/example_slides.pdf"