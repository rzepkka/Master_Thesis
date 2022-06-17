import pandas as pd
import numpy as np
import pickle
import json
import warnings
import os, glob
import subprocess
from PIL import Image
warnings.filterwarnings("ignore")

from mapping_3D import save_subtype_data
from visualize import get_labels, plot_dk_atlas, plot_aseg_atlas

# =================== 1. LOAD PICKLE FILE ========================================================================================================================================================
read_input_file = open('data/EDADS_subtype_timelines_agecorrected_opt.pickle','rb')
load_inputs = pickle.load(read_input_file)
read_input_file.close()

T, S, Sboot = load_inputs

# =================== 2. GET SUBTYPE LABELS (1.) OR SPECIFY YOUR OWN LABELS (2.) =====================================================================================================================================

# 1. 
labels = get_labels(S=S)
# labels = np.array(labels)
# np.array(labels).dump(open('temp_folder/csv/labels.npy', 'wb'))

# 2. 
# labels = ['A', 'B', 'C', 'D']
# labels = np.array(labels)
# np.array(labels).dump(open('temp_folder/csv/labels.npy', 'wb'))
# labels=list(labels)

# =================== 3. LOAD JSON FILES FOR BRAIN REGION MAPPINGS ===============================================================================================

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


# =================== 4. SET UP R ========================================================================================================================
command = 'Rscript'
arg = '--vanilla'
path2script = 'subprocess.R'

def safe_html(labels, command, path2script):
    for i in range(len(labels)):
        i=str(i)
        output = subprocess.run([command, path2script,i])

    print('PROGRESS: All HTML files successfully saved in: /temp_folder')
      # print('CONTINUE running Streamlit APP by calling ">> streamlit run app.py" from command line')


# =================== 6. SAVE 2D ANIMATIONS ========================================================================================================================

def make_gif(frame_folder, subtype, atlas):
    file_list = glob.glob(f'{frame_folder}/*.png')
    file_list.sort()
    frames = [Image.open(image) for image in file_list]
    frame_one = frames[0]
    frame_one.save(f"temp_folder/2D_animations/{atlas}-{subtype}.gif", format="GIF", append_images=frames,
               save_all=True, duration=200, loop=0) 

def generate_animations(T, S, labels, map_dk, map_aseg):

    for subtype in labels:

        # Clear folder
        directory = 'temp_folder/snapshots'
        filelist = glob.glob(os.path.join(directory, "*"))
        for f in filelist:
            os.remove(f)

        video_slider = np.linspace(0,1,50)

        for value in video_slider:
            filename = f"DK-{subtype}-{value}"

            plot_dk_atlas(T = T, S = S, map_dk = map_dk, 
                          subtype = subtype, 
                          slider=value,
                          save=True, 
                          filename=filename)

        make_gif("temp_folder/snapshots", subtype,'DK')

        # Clear folder
        directory = 'temp_folder/snapshots'
        filelist = glob.glob(os.path.join(directory, "*"))
        for f in filelist:
            os.remove(f)

        for value in video_slider:
            filename = f"ASEG-{subtype}-{value}"

            plot_aseg_atlas(T = T, S = S, map_aseg = map_aseg, 
                          subtype = subtype, 
                          slider=value,
                          save=True, 
                          filename=filename)

        make_gif("temp_folder/snapshots", subtype,'ASEG')

        print(f'PROGRESS: animations for {subtype} SUCCESSFULLY saved to folder: temp_folder/2D_animations')



# =================== 5. SET-UP STREAMLIT APP ========================================================================================================================

# RUN MAPPING FUNCTIONS - dk_3D() AND aseg_3D() FOR EACH Subtype AND SAVE DATA IN /temp_folder/csv
save_subtype_data(T, S, map_dk, map_aseg)

# SAVE 3D BRAIN VISUALIZATIONS IN /temp_folder/3D_files
safe_html(labels, command, path2script)

# GENERATE AND SAVE 2D GIF ANIMATIONS IN /temp_folder/2D_animations
generate_animations(T=T, S=S, labels=labels, map_dk=map_dk_2D, map_aseg=map_aseg_2D)

print('CONTINUE running Streamlit APP by calling ">> streamlit run app.py" from command line')










