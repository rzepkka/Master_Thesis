import streamlit as st
import streamlit.components as stc
import streamlit.components.v1 as components

import scipy.stats 
import pandas as pd
from matplotlib import rc
import pickle
# from bokeh.models import Circle, HoverTool, PrintfTickFormatter
# from bokeh.models import LassoSelectTool, WheelZoomTool, PanTool, BoxZoomTool, SaveTool, ResetTool, ColumnDataSource
# import markdown

from streamlit_graphs import event_centers, plot_ggseg, plot_dk_atlas, plot_aseg_atlas, subtypes_pieplot
# from streamlit_echarts import st_echarts

# from mapping_2D import mapping_dk, dk_dict, aseg_dict

from pathlib import Path
import base64
import time
import matplotlib.pyplot as plt

# rpy2
import rpy2
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr, data

import rpy2.robjects as robjects
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri

from IPython.display import HTML


utils = importr('utils')
base = importr("base")
datasets = importr('datasets')
ggseg = importr("ggseg")
ggplot2 = importr("ggplot2")
dplyr = importr("dplyr")
tidyr = importr("tidyr")

ggseg3d = importr('ggseg3d')

st.set_option('deprecation.showPyplotGlobalUse', False)

# Import Pikle File
read_input_file = open('data/ADC_FTLD_subtypes_agecorrected_zscore_final.pickle','rb')
load_inputs = pickle.load(read_input_file)
read_input_file.close()

T, S, X = load_inputs


def main():

    st.set_page_config(layout="wide")
    st.sidebar.header("Adjusting Plots")

    # SELECT PLOT
    plot_type_list = ['Boxplots - event centers', 'Ggseg - 3D', 'Pie Plot']
    chosen_plot_type = st.sidebar.selectbox("Select Plot", plot_type_list)

    # CHOOSE WIDTH AND HEIGHT
    chosen_width = st.sidebar.number_input('Select width of the plot in px',value=1000)
    chosen_height = st.sidebar.number_input('Select height of the plot in px',value=800)

    if chosen_plot_type == 'Boxplots - event centers':


        options = ['Subtype 0','Subtype 1', 'Subtype 2', 'Subtype 3', 'Subtype 4']
        num_subtypes = len(options)
        chosen_subtypes = st.multiselect('Choose subtype to compare:', options)
        color_list = []

        default_color_list = ['#0000ff', '#880000', '#ffa07a', '#04977d', '#fd8ef3']

        for idx, subtype in enumerate(chosen_subtypes):
            subtype_color = st.text_input(f'Choose color for {subtype}', value = f'{default_color_list[idx]}',placeholder='e.g. #000000')
            color_list.append(subtype_color)

        # color_map = {chosen_subtypes[i]: color_list[i] for i in range(len(chosen_subtypes))}

        order_list = chosen_subtypes
        orderBy = st.selectbox("Select subtype to order by",order_list)

        eventCenters = event_centers(T = T,
                                    S = S, 
                                    color_list = color_list,
                                    chosen_subtypes = chosen_subtypes,
                                    subtype_labels = options, 
                                    orderBy = orderBy,
                                    width = chosen_width,
                                    height = chosen_height)

        chosen_2D = orderBy        

        col1, col2 = st.columns([2.2,3])

        if chosen_2D != None:
            slider = st.slider(label = 'Choose regions to display', 
                                help = 'Only display regions with ordering <= to the chosen value',
                                min_value = 0, 
                                max_value = len(T.biomarker_labels), 
                                value=len(T.biomarker_labels), step=1)
        
            with col1:
                ggseg_dk = plot_dk_atlas(T = T, S = S, subtype = chosen_2D, slider = slider)
                st.pyplot(ggseg_dk)

            with col2:     
                ggseg_aseg = plot_aseg_atlas(T = T, S = S, subtype = chosen_2D, slider = slider)       
                st.pyplot(ggseg_aseg)


        st.plotly_chart(eventCenters)

    elif chosen_plot_type == 'Ggseg - 3D':

        
        dk_data = pd.read_csv("data/dk_R_subtype0.csv") 
        with localconverter(robjects.default_converter + pandas2ri.converter):
            dk_data_R = robjects.conversion.py2rpy(dk_data)

        dk_3d = ggseg3d.get_atlas('dk_3d', surface = ["LCBC",'inflated'], hemisphere = ['left','right'])

        colors = robjects.r('''
            c('#6e0101','#ffabab')
            ''')
        st.write(colors)

        plot_dk = ggseg3d.ggseg3d(dk_data_R, atlas = dk_3d, surface = ["LCBC",'inflated'], hemisphere = ['left','right'],
               label = 'region', colour = 'p', palette = colors)


        st.write(type(plot_dk))
        st.write(plot_dk)

        # html_string = plot_dk
        # st.markdown(html_string, unsafe_allow_html=True)

        # print(plot_dk)

        components.iframe('/Users/macos/Documents/GitHub/Master_Thesis/notebook examples/html_outputs/plot_dk.html')



    elif chosen_plot_type == 'Pie Plot':

        options = ['Subtype 1', 'Subtype 2', 'Subtype 3', 'Subtype 4']
        num_subtypes = len(options)

        # DISPLAY PIE CHART
        # DEFAULT NUMBERS FOR PIE CHART
        default_numbers = [567, 98, 511, 324]
        fig, out = subtypes_pieplot(default_numbers, options) # labels, pi0_mean, evn_full, evn,num_subtypes, chosen_width, chosen_height, subtype_color
        st.plotly_chart(fig, use_container_width=True)

    #     # BOX PLOT
    #     chosen_subtypes = st.multiselect('Choose subtype to compare:', options)

    #     color_list = []
    #     default_color_list = ['#9baee4', '#880000', '#ffa07a', '#04977d', '#fd8ef3']
    #     for idx, subtype in enumerate(chosen_subtypes):
    #         subtype_color = st.text_input(f'Choose color for {subtype}', value = f'{default_color_list[idx]}',placeholder='e.g. #000000')
    #         color_list.append(subtype_color)

    #     color_map = {chosen_subtypes[i]: color_list[i] for i in range(len(chosen_subtypes))}
    #     order_list = chosen_subtypes
    #     orderBy = st.selectbox("Select subtype to order by",order_list)

    #     eventCenters = event_centers(T = T,
    #                                 S = S, 
    #                                 color_list = default_color_list,
    #                                 chosen_subtypes = chosen_subtypes,
    #                                 subtype_labels = options, 
    #                                 orderBy = orderBy,
    #                                 width = chosen_width,
    #                                 height = chosen_height)
    #     st.plotly_chart(eventCenters)

main()









