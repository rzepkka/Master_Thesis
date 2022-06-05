import pandas as pd
import numpy as np
import pickle
import json
import warnings
import os, glob
import subprocess
warnings.filterwarnings("ignore")

from mapping_3D import dk_3D, aseg_3D, save_subtype_data

def get_labels(S):
      unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
      subtype_labels = []
      for i in range(len(unique_subtypes)):
          subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))        
      return subtype_labels

# =================== 1. LOAD PICKLE FILE ========================================================================================================================================================
read_input_file = open('data/EDADS_subtype_timelines_agecorrected_opt.pickle','rb')
load_inputs = pickle.load(read_input_file)
read_input_file.close()

T, S, Sboot = load_inputs

# =================== 2. GET SUBTYPE LABELS=====================================================================================================================================

labels = get_labels(S=S)

# =================== 3. LOAD JSON FILES FOR BRAIN REGION MAPPINGS ===============================================================================================
f = open('data/DK_3D_combined.json')
map_dk = json.load(f)
f.close()

f = open('data/ASEG_3D_combined.json')
map_aseg = json.load(f)
f.close()

# CLEAN TEMP FOLDER
# directory = 'temp_folder'
# filelist = glob.glob(os.path.join(directory, "*"))
# for f in filelist:
#     os.remove(f)

# =================== 4. RUN MAPPING FUNCTIONS - dk_3D() AND aseg_3D() FOR EACH Subtype AND SAVE DATA IN temp_folder
save_subtype_data(T, S, map_dk, map_aseg)


# =================== 5. SET UP R ========================================================================================================================
command = 'Rscript'
arg = '--vanilla'
path2script = 'subprocess.R'

def safe_html(labels, command, path2script):
    for i in range(len(labels)):
        i=str(i)
        output = subprocess.run([command, path2script,i])

    print('PROGRESS: All HTML files successfully saved in: /temp_folder')
    print('CONTINUE running Streamlit APP by calling ">> streamlit run app.py" from command line')

safe_html(labels, command, path2script)















