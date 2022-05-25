import pandas as pd
import numpy as np
import pickle
import json
import warnings
warnings.filterwarnings("ignore")

from mapping_3D import dk_3D, aseg_3D, save_subtype_data

# 1. LOAD PICKLE FILE
read_input_file = open('data/EDADS_subtype_timelines_agecorrected_opt.pickle','rb')
load_inputs = pickle.load(read_input_file)
read_input_file.close()

T, S, Sboot = load_inputs

# 2. Load JSON files
f = open('data/DK_3D_combined.json')
map_dk = json.load(f)
f.close()

f = open('data/ASEG_combined.json')
map_aseg = json.load(f)
f.close()

# 3. Run mapping functions - dk_3D() and aseg_3D() for each Subtype and save data in temp folder


save_subtype_data(T, S, map_dk, map_aseg)


