import numpy as np
import scipy.stats 
import pandas as pd
import os, glob
from PIL import Image



def save_subtype_data(T, S, map_dk, map_aseg):  
    """
    Saves all subtypes data in separate .csv files for DK and ASEG atlas
    :param T: Timeline object
    :param S: subtyping dictionary, subtypes for each patient individually
    :param map_dk: dictionary from loaded JSON file, key: values --> T.biomarker_labels: list(DK-labels)
    :param map_aseg: dictionary from loaded JSON file, key: values --> T.biomarker_labels: list(ASEG-labels)
    :return: info if data was saved
    """
    
    # Get subtype labels
    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    subtype_labels = []
    for i in range(len(unique_subtypes)):
        subtype_labels.append('Subtype '+str(int(unique_subtypes[i])))  
        
    # Save each file as csv
    for label in subtype_labels:
        
        dk = dk_3D(T, S, 
                   mapped_dict = map_dk, 
                   subtype = label)
        
        aseg = aseg_3D(T, S, 
                       mapped_dict = map_aseg, 
                       subtype = label) 

        dk.to_csv(f'temp_folder/csv/dk_R_{label}.csv', index = False)
        aseg.to_csv(f'temp_folder/csv/aseg_R_{label}.csv', index = False)    
        
    print('PROGRESS: All suptype files saved in: /temp_folder/csv')


def safe_html(labels, command, path2script):
    for i in range(len(labels)):
        i=str(i)
        output = subprocess.run([command, path2script,i])

    print('PROGRESS: All HTML files successfully saved in: /temp_folder')


# ====================================== DK-ATLAS ==================================================================================================================

def dk_3D(T,S, mapped_dict, subtype_labels = None, subtype = None):
    
    """
    Creates a dictionary, which can be used as input to ggseg3d() function
    :param T: Timeline object
    :param S: subtyping dictionary, subtypes for each patient individually
    :param mapped_dict: dictionary from loaded JSON file, key: values --> T.biomarker_labels: list(DK-labels)
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
    labels = list(map(lambda x: x.lower(), T.biomarker_labels))
    labels_dict = {num: label.lower() for num, label in enumerate(labels)}
    
    order = T.sequence_model['ordering'][subtype_labels[subtype]]
    
    labels_ordered = []
    for o in order:
        labels_ordered.append(labels_dict[o])       
    
    dic = dict(zip(labels_ordered, T.sequence_model['event_centers'][subtype_labels[subtype]]))
                    
    # flat lost of dict values (single list of DK-labels)
    dk_flat = [x for v in mapped_dict.values() for x in v]
        
    hemi = []
    for idx, region in enumerate(dk_flat):
        if 'left' in region:
            hemi.append('left')
            dk_flat[idx]=dk_flat[idx].replace(' left','')
        elif 'right' in region:
            hemi.append('right')
            dk_flat[idx]=dk_flat[idx].replace(' right','')
        else:
            hemi.append('subcort')
              
    #Match T.biomarker_labels to DK labels
    list_plot = list()
    for key in mapped_dict.keys():
        for item in mapped_dict[key]:
            list_plot.append(dic[key])
    
    dic_dk = {'region': dk_flat, 'hemi':hemi, 'p': list_plot}
    df = pd.DataFrame(dic_dk)
    
    return df


# ====================================== ASEG atlas ==================================================================================================================

def aseg_3D(T, S, mapped_dict, subtype_labels = None, subtype = None):
    
    """
    Creates a dictionary, which can be used as input to ggseg.plot_dk() function
    :param T: Timeline object
    :param S: subtyping dictionary, subtypes for each patient individually
    :param mapped_dict: dictionary from loaded JSON file, key: values --> T.biomarker_labels: list(DK-labels)
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
    dic_aseg = {'region': aseg_flat, 'p':list_plot} 
    
    df = pd.DataFrame(dic_aseg)
        
    return df