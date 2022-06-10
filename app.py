import streamlit as st
import streamlit.components as stc
import streamlit.components.v1 as components

import scipy.stats 
import pandas as pd
import numpy as np
import pickle
import numpy as np
import os
import webbrowser
from collections import Counter
import collections
import json
import re
import glob
from PIL import Image


import matplotlib.pyplot as plt
from matplotlib import rc

import plotly.graph_objs as go
from ipywidgets import Output

from visualize import event_centers, plot_dk_atlas, plot_aseg_atlas, patient_staging, staging_boxes

from visualize import atypicality, atypicality_boxes, staging_scatterplot, piechart_multiple, custom_dk, custom_aseg, plot_ggseg

from visualize import subtype_probabilities, individual_staging, biomarker_distribution


st.set_option('deprecation.showPyplotGlobalUse', False)

# ===================== LOAD FILES ==============================================================================================================================
# LOAD PICKLE FILE
read_input_file = open('data/EDADS_subtype_timelines_agecorrected_opt.pickle','rb')
load_inputs = pickle.load(read_input_file)
read_input_file.close()

T, S, Sboot = load_inputs

# LOAD DIAGNOSIS VARIABLE
diagnosis = np.load('data/diagnosis.npy', allow_pickle=True)

# LOAD JSON FILES FOR MAPPING BRAIN REGIONS
f = open('data/DK_2D_combined.json')
map_dk = json.load(f)
f.close()

f = open('data/ASEG_2D_combined.json')
map_aseg = json.load(f)
f.close()

# Get labels for options in select boxes
def get_labels(S):
      unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
      subtype_labels = []
      for i in range(len(unique_subtypes)):
          subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))        
      return subtype_labels

# Create 2D gifs
def make_gif(frame_folder, subtype, atlas):
    file_list = glob.glob(f'{frame_folder}/*.png')
    file_list.sort()
    frames = [Image.open(image) for image in file_list]
    frame_one = frames[0]
    frame_one.save(f"animations/{atlas}-{subtype}.gif", format="GIF", append_images=frames,
               save_all=True, duration=200, loop=0) 

labels = get_labels(S=S)

# Setup the APP - saving all files 

# ===================== APP ==============================================================================================================================
def main():

    st.set_page_config(layout="wide")
    st.sidebar.title("Adjusting Plots")

    # Connect .css file for styling the app components
    def local_css(file_name):
        with open(file_name) as f:
            st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
    local_css("style.css")

    # SELECT PLOT
    plot_type_list = ['Disease timeline','Individual', 'Project overview']
    chosen_plot_type = st.sidebar.selectbox("Select Plot", plot_type_list)


    if chosen_plot_type == 'Disease timeline':

        st.header('Disease progression timeline')

        # CHOOSE WIDTH AND HEIGHT
        chosen_width = st.sidebar.number_input('Select width of the plot in px',value=1000, max_value=1130)
        chosen_height = st.sidebar.number_input('Select height of the plot in px',value=800)

        col_piechart, col_piechart_select = st.columns([5.2,3])

        with col_piechart:
            diagnosis_labels = list(set(diagnosis))
            diagnosis_labels.remove('CN')
            choose_pieplot = st.multiselect('Diagnoses included in the piechart:', diagnosis_labels, default=diagnosis_labels)
            plot_piechart = piechart_multiple(S=S,
                                            diagnosis=diagnosis,
                                            chosen_subtypes=choose_pieplot)
            st.plotly_chart(plot_piechart)

        # SELECTION
        num_subtypes = len(labels)
 
        subtype_list = []

        col_select, col_select_blank = st.columns([5.2,3])
        with col_select:
            subtype_visualize = st.selectbox('Select a subtype to visualize:',labels)       
        subtype_list.append(subtype_visualize) 

        # ======================= 2D ===============================================================================================================
        
        col_cortical, col_subcortical, col_button = st.columns([2,3.2,3])

        if subtype_visualize != None:

            # Clear folder
            directory = 'video'
            filelist = glob.glob(os.path.join(directory, "*"))
            for f in filelist:
                os.remove(f)

            col_slider, col_slider_blank = st.columns([5.2,3])
            with col_slider:
                slider = st.slider(label = 'Choose regions to display', 
                                    help = 'Only display regions with ordering <= to the chosen value',
                                    min_value = 0.0, 
                                    max_value = 1.0, 
                                    value=1.0, step=0.01)
        
            with col_cortical:
                ggseg_dk = plot_dk_atlas(T = T, 
                                        S = S,
                                        map_dk=map_dk,
                                        subtype = subtype_visualize, 
                                        slider = slider)

                st.pyplot(ggseg_dk)

            with col_subcortical:     
                ggseg_aseg = plot_aseg_atlas(T = T, 
                                        S = S, 
                                        map_aseg=map_aseg,
                                        subtype = subtype_visualize, 
                                        slider = slider)       
                st.pyplot(ggseg_aseg)

            # BUTTONS
            html_file = subtype_visualize
            # html_file = subtype_visualize.replace(" ","_")

            color_list = []
            default_color_list = ['#1f77b4', '#ff7f0e', '#2ca02c','#d62728', '#9467bd'] 

            with col_button:
                if st.button('Open 3D visualization in a new tab'):
                    try:
                        filename = "file://"+os.getcwd()+ f"/temp_folder/{html_file}.html"
                        webbrowser.open(filename)

                    except FileNotFoundError:
                        st.sidebar.error('File not found')
                        st.sidebar.error('Please run app_setup file before trying to download 3D visualisations')

# ======================= ANIMATIONS ===============================================================================================================

                # DK animation
                if st.button('Generate 2D animations'):

                    # Clear folder
                    directory = 'video'
                    filelist = glob.glob(os.path.join(directory, "*"))
                    for f in filelist:
                        os.remove(f)

                    video_slider = np.linspace(0,1,50)

                    for value in video_slider:
                        filename = f"{subtype_visualize}-{value}"
                        
                        plot_dk_atlas(T = T, S = S, map_dk = map_dk, 
                                      subtype = subtype_visualize, 
                                      slider=value,
                                      save=True, 
                                      filename=filename)

                    make_gif("video", subtype_visualize,'DK')

                    # Clear folder
                    directory = 'video'
                    filelist = glob.glob(os.path.join(directory, "*"))
                    for f in filelist:
                        os.remove(f)

                    for value in video_slider:
                        filename = f"ASEG-{subtype_visualize}-{value}"
                        
                        plot_aseg_atlas(T = T, S = S, map_aseg = map_aseg, 
                                      subtype = subtype_visualize, 
                                      slider=value,
                                      save=True, 
                                      filename=filename)  # , save=True

                    make_gif("video", subtype_visualize,'ASEG')


                    with col_cortical:
                        st.image(f"animations/DK-{subtype_visualize}.gif")
                    with col_subcortical:
                        st.image(f"animations/ASEG-{subtype_visualize}.gif")


        # ======================= EVENT CENTERS ===============================================================================================================

        # list for additional subtypes to compare
        options_compare = []
        for label in labels:
            if label == subtype_visualize:
                pass
            else: 
                options_compare.append(label)

        col_event_centers, col_event_centers_options = st.columns([3,1])

        with col_event_centers_options:
            st.subheader('Style event centers graph')
            chosen_subtypes = st.multiselect('Select additional subtypes to compare:', options_compare)

            title_font = st.number_input('Title size:',value=34)
            title_axes = st.number_input('Axis labels size:',value=18)
            title_ticks = st.number_input('Ticks size:',value=14)
            title_legend = st.number_input('Legend size:',value=22)

            font_list = [title_font, title_axes, title_ticks, title_legend]

            for subtype in chosen_subtypes:
                subtype_list.append(subtype)

            for idx, subtype in enumerate(subtype_list):
                subtype_color = st.text_input(f'Select color for {subtype}', value = f'{default_color_list[idx]}',placeholder='e.g. #000000')

                match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', subtype_color)

                if match:
                    color_list.append(subtype_color)
                else:
                    color_list.append(default_color_list[idx])
                    st.error('Please specify a valid hex volor value, e.g. #000000.')


        eventCenters = event_centers(T = T,
                                    S = S, 
                                    color_list = color_list,
                                    chosen_subtypes = subtype_list,
                                    subtype_labels = labels, 
                                    orderBy = subtype_visualize,
                                    width = chosen_width,
                                    height = chosen_height,
                                    slider = slider,
                                    fontsize = font_list)

    
        with col_event_centers:
            st.plotly_chart(eventCenters)

        # ADD DIVIDER
        st.markdown('---')

        # ======================= PATIENT STAGING ===============================================================================================================

        col_staging, col_staging_options = st.columns([3,1])

        diagnosis_labels = list(set(diagnosis))
        diagnosis_labels.sort()

        color_list = ['#e41a1c','#377eb8','#0000ff','#00D612','#1f77b4', '#ff7f0e', '#2ca02c','#d62728', '#9467bd']
        color_diagnosis =[]

        with col_staging:
            chosen_plot = st.selectbox('Select plot:',['Patient Staging', 'Scatterplot'])

        #  ================== PATIENT STAGING =================================================================================================================

        with col_staging_options:
            st.subheader('Style the graphs') 
            opacity = st.number_input('Select opacity value:', value=0.8)
            title_font = st.number_input('Title font:',value=34)
            title_axes = st.number_input('Axis labels font:',value=18)
            title_ticks = st.number_input('Ticks size font:',value=14)
            title_legend = st.number_input('Legend font:',value=22)
            font_list = [title_font, title_axes, title_ticks, title_legend]


        if chosen_plot =='Patient Staging':

            with col_staging_options:

                num_bins = st.number_input('Select the number of bins', value = 10)
                bin_width = st.number_input('Select bin width:', value = 0.02)

                for idx, label in enumerate(diagnosis_labels):
                    if label != 'CN':
                        color = st.text_input(f'Select color for {label}', value = f'{color_list[idx]}',placeholder='e.g. #000000')


                        match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', color)

                        if match:
                            color_diagnosis.append(color)
                        else:
                            color_diagnosis.append(color_list[idx])
                            st.error('Please specify a valid hex volor value, e.g. #000000.')
                        

            with col_staging:
                barmode = st.radio('Select barmode:', ['group','stack'])

                # BARPLOT
                plot_staging = patient_staging(S=S,
                                        diagnosis=diagnosis, 
                                        color_list = color_diagnosis,
                                        num_bins=num_bins, 
                                        bin_width=bin_width,
                                        width = chosen_width,
                                        height = chosen_height,
                                        fontsize=font_list,
                                        opacity=opacity,
                                        barmode=barmode)

                # BOX
                box_staging = staging_boxes(S=S,
                                        diagnosis=diagnosis,
                                        color_list=color_diagnosis,
                                        width=chosen_width,
                                        fontsize=font_list)


                st.plotly_chart(plot_staging)
                st.plotly_chart(box_staging)  

      
        if chosen_plot =='Scatterplot':

            color_list = ['#1f77b4', '#ff7f0e', '#2ca02c','#d62728', '#9467bd']  

            with col_staging_options:

                for idx, label in enumerate(labels):
                    color = st.text_input(f'Color for {label}', value = f'{color_list[idx]}',placeholder='e.g. #000000')
                    color_diagnosis.append(color)


            with col_staging:
                # diagnosis_labels = list(set(diagnosis))
                # diagnosis_labels.remove('CN')
                choose_scatterplot = st.multiselect('Diagnoses included in the scatterplot:', labels, default=labels)

                plot_scatter = staging_scatterplot(S=S,
                                                diagnosis=diagnosis,
                                                chosen_subtypes=choose_scatterplot,
                                                color_list=color_diagnosis,
                                                width = chosen_width,
                                                height = chosen_height)

                st.plotly_chart(plot_scatter)

        # ADD DIVIDER
        st.markdown('---')

    elif chosen_plot_type == 'Individual':

        data = pd.read_csv("data/EDADS_data.csv")

        st.header('Patient-specific information')

        # CHOOSE WIDTH AND HEIGHT
        chosen_width = st.sidebar.number_input('Select width of the plot in px',value=900, max_value=1130)
        chosen_height = st.sidebar.number_input('Select height of the plot in px',value=600)

        patient_id = int(st.text_input(f"Select patient's ID", value = 0,placeholder='id...'))

        col_probabilities, col_probabilities_options = st.columns([3,1])


        with col_probabilities_options:
            st.subheader('Style the graphs') 
            title_font = st.number_input('Title font:',value=34)
            title_axes = st.number_input('Axis labels font:',value=18)
            title_ticks = st.number_input('Ticks size font:',value=14)
            title_bars = st.number_input('Annotation font:',value=22)
            font_list = [title_font, title_axes, title_ticks, title_bars]


            color = st.text_input(f'Select color', value = '#1f77b4',placeholder='e.g. #000000')

            match = re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', color)

            if match:
                pass
            else:
                color = '#636EFA'
                st.error('Please specify a valid hex volor value, e.g. #000000.')

        with col_probabilities:

            # st.subheader(f"Patients diagnosis: {data['Diagnosis'][data['PTID']==patient_id][0]}")
            # st.subheader(f"Patients prediction: {prediction}")

           

            if patient_id in list(data['PTID']):

                plot_probabilities, prediction = subtype_probabilities(info=data,
                                                                S=S,
                                                                patient_id=patient_id,
                                                                fontlist=font_list,
                                                                color=color,
                                                                width = chosen_width,
                                                                height = chosen_height
                                                                )

                d = np.array(data['Diagnosis'][data['PTID']==patient_id])[0]
                st.subheader(f"Patients diagnosis: {d}")
                st.subheader(f"Patients prediction: {prediction}")

                box_individual = individual_staging(data=data,
                                                Sboot=Sboot,
                                                patient_id=patient_id,
                                                color=color,
                                                fontsize=font_list[1],
                                                width=chosen_width
                                                )

                st.plotly_chart(plot_probabilities)
                st.plotly_chart(box_individual)

                if prediction not in  ['Outlier', 'No prediction']:
                    subtype = int(prediction[-1])
                else:
                    subtype = 0

                plot_distribution = biomarker_distribution(data=data,
                                            T=T,
                                        subtype=subtype,
                                        patient_id=patient_id)

                st.plotly_chart(plot_distribution)

            else: 
                st.info("Provide an existing patient's ID")

        #     if prediction not in  ['Outlier', 'No prediction']:
        #         subtype = int(prediction[-1])
        #     else:
        #         subtype = 0

        #     plot_distribution = biomarker_distribution(data=data,
        #                                 T=T,
        #                                 subtype=subtype,
        #                                 patient_id=patient_id)

        # st.plotly_chart(plot_distribution)
        

    else:
        st.subheader('Project presentation will be here')

main()









