import numpy as np
import scipy.stats 
import pandas as pd
import os, glob
from PIL import Image

# from visualize import plot_dk_atlas, plot_aseg_atlas

# ====================================== NEW DATA ==================================================================================================================
# FUNCTIONS FOR T.BIOMARKER_LABELS WITH BRAIN REGIONS CONSIDERED SEPARATELY FOR LEFT AND RIGHT HEMISPHERE

def dk_dict(T,S, mapped_dict, subtype_labels = None, subtype = None):
    
    """
    Creates a dictionary, which can be used as input to ggseg.plot_dk() function
    :param T: dataframe from dk_dataframe() function
    :param S: chosen subtype
    :param mapped_dict: a dictionary with key: values --> T.biomarker_labels: list(DK-labels)
    :param subtype: name or index of the subtype from subtype_lables (optional, choses first available subtype as default)  
    :param subtype_labels: a list with names of the subtypes (optional)
    :return: dictionary with scores for each DK region for chosen subtype
    """
    
    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    if subtype_labels is None:
        subtype_labels = {f'Subtype {i}': i for i in range(len(unique_subtypes))}
        if subtype is None:
            subtype = next(iter(subtype_labels))
    elif subtype is None:
        subtype = subtype_labels[0]  
        
        
    # clean names from capital letters
    labels = list(map(lambda x: x.lower(), T.biomarker_labels))
    labels_dict = {num: label.lower() for num, label in enumerate(labels)}
    
    order = T.sequence_model['ordering'][subtype_labels[subtype]]
    
    labels_ordered = []
    for o in order:
        labels_ordered.append(labels_dict[o])   
    
    # All regions
    dic = dict(zip(labels_ordered, T.sequence_model['event_centers'][subtype_labels[subtype]]))
                    
    # flat list of dict values (single list of DK-labels)
    dk_flat = [x for v in mapped_dict.values() for x in v]
                
    #Match T.biomarker_labels to DK labels
    list_plot = list()
    for key in mapped_dict.keys():
            for item in mapped_dict[key]:
                list_plot.append(dic[key])
                    
    # Dict for dk-label: T.label value
    dic_dk = dict(zip(dk_flat, list_plot))
    
    return dic_dk


def aseg_dict(T, S, mapped_dict, subtype_labels = None, subtype = None, hemi = False):
    
    """
    Creates a dictionary, which can be used as input to ggseg.plot_dk() function
    :param T: Timeline object
    :param S: subtyping dictionary, subtypes for each patient individually
    :param mapped_dict: a dictionary with key: values --> T.biomarker_labels: list(DK-labels)
    :param subtype_labels: a list with names of the subtypes (optional)
    :param subtype: name or index of the subtype from subtype_lables (optional, choses first available subtype as default)  
    :return: dictionary with scores for each DK region for chosen subtype
    """

    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    if subtype_labels is None:
        subtype_labels = {f'Subtype {i}': i for i in range(len(unique_subtypes))}
        if subtype is None:
            subtype = next(iter(subtype_labels))
    elif subtype is None:
        subtype = subtype_labels[0]
    
    # clean names from capital letters
    labels = T.biomarker_labels
    labels_dict = {num: label.lower() for num, label in enumerate(labels)}
    order = T.sequence_model['ordering'][subtype_labels[subtype]]
    
    labels_ordered = []
    for o in order:
        labels_ordered.append(labels_dict[o])    
     
    # Dictionary with all labels
    dic = dict(zip(labels_ordered, T.sequence_model['event_centers'][subtype_labels[subtype]]))
    
    # flat list of dict values (single list of DK-labels)
    aseg_flat = [x for v in mapped_dict.values() for x in v]
            
    #Match T.biomarker_labels to DK labels
    list_plot = list()
    for key in mapped_dict.keys():
        for item in mapped_dict[key]:
            list_plot.append(dic[key])
            
    # Dict for dk-label: T.label value
    dic_aseg = dict(zip(aseg_flat, list_plot))    
        
    return dic_aseg






