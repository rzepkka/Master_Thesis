import streamlit as st
import streamlit.components as stc
import streamlit.components.v1 as components

import scipy.stats 
import pandas as pd
import numpy as np
import numpy as np
import os
import webbrowser
import re
import glob
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import rc
import base64
import plotly.graph_objs as go
from ipywidgets import Output

from visualize import event_centers, plot_dk_atlas, plot_aseg_atlas, patient_staging, staging_boxes
from visualize import piechart_multiple, custom_dk, custom_aseg
from visualize import subtype_probabilities, individual_staging, biomarker_distribution
from visualize import get_prediction, make_gif, displayPDF

from input_data import data
from input_data import T, S, Sboot, diagnosis
from input_data import map_dk_2D as map_dk
from input_data import map_aseg_2D as map_aseg
from input_data import subtype_labels
from input_data import disease


st.set_option('deprecation.showPyplotGlobalUse', False)


# ===================== STREAMLIT APP ==============================================================================================================================
def main(disease, subtype_labels):

    st.set_page_config(layout="wide")
    st.sidebar.title("Menu")

    # Specify subtype labels
    # subtype_labels = ['Subcortical subtype', 'Frontal subtype', 'Parietal subtype','Typical subtype']

    # Connect .css file for styling the app components
    def local_css(file_name):
        with open(file_name) as f:
            st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
    local_css("style.css")

    # SELECT PLOT
    plot_type_list = [f'Disease timeline for {disease}','Patient-specific information','Load a report']
    chosen_plot_type = st.sidebar.radio("", plot_type_list)

    if chosen_plot_type == f'Disease timeline for {disease}':

        st.header(f'Disease progression timeline for {disease}')

        # CHOOSE WIDTH AND HEIGHT
        chosen_width = st.sidebar.number_input('Select width of the plot in px',value=1000, max_value=1130)
        chosen_height = st.sidebar.number_input('Select height of the plot in px',value=800)

        col_piechart, col_piechart_select = st.columns([5.2,3])

        with col_piechart:
            diagnosis_labels = list(set(diagnosis))
            diagnosis_labels.remove('CN')
            choose_pieplot = st.multiselect('Diagnoses included in the piechart:', diagnosis_labels, default=diagnosis_labels,
                                            help ='Pie chart present the distribution of patients with respect to different diagnoses and subtypes of the disease.')
            labels = subtype_labels.copy()

            plot_piechart = piechart_multiple(S=S,
                                            diagnosis=diagnosis,
                                            chosen_subtypes=choose_pieplot,
                                            subtype_labels=labels
                                            )
            st.plotly_chart(plot_piechart)

        # SELECTION
        # subtype_labels = ['Subcortical subtype', 'Frontal subtype', 'Parietal subtype','Typical subtype']
        num_subtypes = len(subtype_labels)
        subtype_list = []
        col_select, col_select_blank = st.columns([5.2,3])
        with col_select:
            subtype_visualize = st.selectbox('Select a subtype to visualize:',subtype_labels)       
        subtype_list.append(subtype_visualize) 

# ======================= 2D STATIC ===============================================================================================================
        
        chosen_2D = st.radio('2D display', ['Static','Animation'], 
                            help = '2D brain visualisations of the sequence in which biomarkers for a particular disease subtype become abnormal, for both cortical and subcortical regions of the brain.')
        col_cortical, col_subcortical, col_button = st.columns([2,3.2,3])

        if subtype_visualize != None:
            col_slider, col_slider_blank = st.columns([5.2,3])
            with col_slider:
                slider = st.slider(label = 'Choose regions to display', 
                                    help = 'Only display regions with ordering <= to the chosen value',
                                    min_value = 0.0, 
                                    max_value = 1.0, 
                                    value=1.0, step=0.01)

            if chosen_2D == 'Static':
                subtype_labels = ['Subcortical subtype', 'Frontal subtype', 'Parietal subtype','Typical subtype']
        
                with col_cortical:
                    ggseg_dk = plot_dk_atlas(T = T, 
                                            S = S,
                                            map_dk=map_dk,
                                            subtype = subtype_visualize, 
                                            slider = slider,
                                            subtype_labels=subtype_labels)

                    st.pyplot(ggseg_dk)

                with col_subcortical:     
                    ggseg_aseg = plot_aseg_atlas(T = T, 
                                            S = S, 
                                            map_aseg=map_aseg,
                                            subtype = subtype_visualize, 
                                            slider = slider,
                                            subtype_labels=subtype_labels)       
                    st.pyplot(ggseg_aseg)

# ======================= 2D ANIMATIONS ===============================================================================================================

            elif chosen_2D == 'Animation':

                with col_cortical:
                    st.image(f"temp_folder/2D_animations/DK-{subtype_visualize}.gif")
                with col_subcortical:
                    st.image(f"temp_folder/2D_animations/ASEG-{subtype_visualize}.gif")

            # BUTTONS
            html_file = subtype_visualize

            color_list = []
            default_color_list = ['#1f77b4', '#ff7f0e', '#2ca02c','#d62728', '#9467bd'] 

            with col_button:
                if st.button('Open 3D visualization in a new tab'):
                    try:
                        filename = "file://"+os.getcwd()+ f"/temp_folder/3D_files/{html_file}.html"
                        webbrowser.open(filename)

                    except FileNotFoundError:
                        st.sidebar.error('File not found')
                        st.sidebar.error('Please run app_setup file before trying to download 3D visualisations')

# ======================= EVENT CENTERS ===============================================================================================================

        # list for additional subtypes to compare
        options_compare = []
        for label in subtype_labels:
            if label == subtype_visualize:
                pass
            else: 
                options_compare.append(label)

        col_event_centers, col_event_centers_options = st.columns([3,1])

        with col_event_centers_options:
            st.subheader('Style event centers graph')
            chosen_subtypes = st.multiselect('Select additional subtypes to compare:', options_compare,
                                                help='Boxplot of event centers for each disease subtype, which illustrates the estimated order of the occurring abnormalities, measured by a 100 repetitions of bootstrapping.')

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
                                    subtype_labels = subtype_labels, 
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

        with col_staging_options:
            st.subheader('Style the graphs') 
            opacity = st.number_input('Select opacity value:', value=0.8)
            title_font = st.number_input('Title font:',value=34)
            title_axes = st.number_input('Axis labels font:',value=18)
            title_ticks = st.number_input('Ticks size font:',value=14)
            title_legend = st.number_input('Legend font:',value=22)
            font_list = [title_font, title_axes, title_ticks, title_legend]

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
            barmode = st.radio('Select barmode:', ['group','stack'],
                                help='Patient staging plots that present the distribution of patients within each diagnostic group, with respect to the predicted disease stage.')

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

        # ADD DIVIDER
        st.markdown('---')

    elif chosen_plot_type == 'Patient-specific information':

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
        

            if patient_id in list(data['PTID']):

                plot_probabilities, prediction = subtype_probabilities(info=data,
                                                                S=S,
                                                                patient_id=patient_id,
                                                                fontlist=font_list,
                                                                color=color,
                                                                width = chosen_width,
                                                                height = chosen_height,
                                                                subtype_labels=subtype_labels
                                                                )

                prediction = get_prediction(data=data,
                                             S=S,
                                              patient_id=patient_id,
                                             subtype_labels=subtype_labels)

                d = np.array(data['Diagnosis'][data['PTID']==patient_id])[0]
                st.subheader(f"Patient's diagnosis: {d}")
                st.subheader(f"Predicted subtype: {prediction}")

                color_list = ['#e41a1c','#377eb8','#4daf4a','#FFA500']

                box_individual = individual_staging(data=data,
                                                S=S,
                                                Sboot=Sboot,
                                                patient_id=patient_id,
                                                diagnosis=diagnosis,
                                                color_list=color_list,
                                                height=600,
                                                width=900)

                st.plotly_chart(plot_probabilities)
                st.plotly_chart(box_individual)

                unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
                s = {subtype_labels[i]: i for i in range(len(unique_subtypes))}
                subtype = s[prediction]

                st.markdown('---')

                plot_distribution = biomarker_distribution(data=data,
                                            T=T,
                                        subtype=subtype,
                                        patient_id=patient_id)

                st.plotly_chart(plot_distribution)

            else: 
                st.info("Provide an existing patient's ID")
  
    else:
        st.subheader('Presentation slides')

        path_to_pdf = st.text_input(f'Input path to PDF file and press Enter to apply: ', value = f'',placeholder='...')

        if path_to_pdf != "":

            displayPDF(path_to_pdf)


# ========== RUN THE APP ==================================================================================================================================

main(disease, subtype_labels)









