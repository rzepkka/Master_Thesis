import numpy as np
import scipy.stats 
import pandas as pd

# ====================================== NEW DATA ==================================================================================================================

def mapping_dk(hemi = False):
    # if hemi = True, regions are mapped separately to left and right hemisphere
    
    # LEFT
    list_cortical_left = ['Temporal_lobe_left','Superior_frontal_gyrus_left',
                          'Middle_frontal_gyrus_left','Inferior_frontal_gyrus_left', 
                          'Gyrus_rectus_left','Orbitofrontal_gyri_left','Precentral_gyrus_left',
                          'Postcentral_gyrus_left','Superior_parietal_gyrus_left', 
                          'Inferolateral_remainder_of_parietal_lobe_left',
                          'Lateral_remainder_of_occipital_lobe_left','Lingual_gyrus_left', 
                          'Insula_left','Gyrus_cinguli_anterior_part_left','Gyrus_cinguli_posterior_part_left',
                          'Parahippocampal_and_ambient_gyri_left']

    org_cortical_mapping_left = [['bankssts_left','transversetemporal_left',
                                      'superiortemporal_left','temporalpole_left','entorhinal_left',
                                      'middletemporal_left','inferiortemporal_left','fusiform_left'], 
                                     ['superiorfrontal_left','frontalpole_left'], 
                                     ['caudalmiddlefrontal_left','rostralmiddlefrontal_left'], 
                                     ['parsopercularis_left','parsorbitalis_left','parstriangularis_left'], 
                                     ['medialorbitofrontal_left'], ['lateralorbitofrontal_left'], 
                                     ['precentral_left','paracentral_left'], ['postcentral_left'], 
                                     ['superiorparietal_left','precuneus_left','cuneus_left','pericalcarine_left'], 
                                 ['inferiorparietal_left','supramarginal_left'], 
                                     ['lateraloccipital_left'], 
                                     ['lingual_left'], ['insula_left'], ['caudalanteriorcingulate_left','rostralanteriorcingulate_leftvolume'], 
                                     ['posteriorcingulate_left','isthmuscingulate_left'], 
                                     ['parahippocampal_left']]
    
    # RIGHT
    list_cortical_right = ['Temporal_lobe_right','Superior_frontal_gyrus_right','Middle_frontal_gyrus_right',
                           'Inferior_frontal_gyrus_right', 'Gyrus_rectus_right','Orbitofrontal_gyri_right',
                           'Precentral_gyrus_right','Postcentral_gyrus_right','Superior_parietal_gyrus_right', 
                           'Inferolateral_remainder_of_parietal_lobe_right','Lateral_remainder_of_occipital_lobe_right',
                           'Lingual_gyrus_right', 'Insula_right','Gyrus_cinguli_anterior_part_right',
                           'Gyrus_cinguli_posterior_part_right','Parahippocampal_and_ambient_gyri_right']
    
    org_cortical_mapping_right = [['bankssts_right','transversetemporal_right',
                                  'superiortemporal_right','temporalpole_right','entorhinal_right',
                                  'middletemporal_right','inferiortemporal_right','fusiform_right'], 
                                 ['superiorfrontal_right','frontalpole_right'], 
                                 ['caudalmiddlefrontal_right','rostralmiddlefrontal_right'], 
                                 ['parsopercularis_right','parsorbitalis_right','parstriangularis_right'], 
                                 ['medialorbitofrontal_right'], ['lateralorbitofrontal_right'], 
                                 ['precentral_right','paracentral_right'], ['postcentral_right'], 
                                 ['superiorparietal_right','precuneus_right','cuneus_right','pericalcarine_right'], ['inferiorparietal_right','supramarginal_right'], 
                                 ['lateraloccipital_right'], 
                                 ['lingual_right'], ['insula_right'], ['caudalanteriorcingulate_right','rostralanteriorcingulate_right'], 
                                 ['posteriorcingulate_right','isthmuscingulate_right'], 
                                 ['parahippocampal_right']]
    
    # BOTH HEMISPHERES
    list_cortical = ['Temporal_lobe','Superior_frontal_gyrus','Middle_frontal_gyrus','Inferior_frontal_gyrus', 
                     'Gyrus_rectus','Orbitofrontal_gyri','Precentral_gyrus','Postcentral_gyrus','Superior_parietal_gyrus', 
                     'Inferolateral_remainder_of_parietal_lobe','Lateral_remainder_of_occipital_lobe',
                     'Lingual_gyrus', 'Insula','Gyrus_cinguli_anterior_part','Gyrus_cinguli_posterior_part',
                     'Parahippocampal_and_ambient_gyri']

    org_cortical_mapping = [['bankssts_left','transversetemporal_left',
                              'superiortemporal_left','temporalpole_left','entorhinal_left',
                              'middletemporal_left','inferiortemporal_left','fusiform_left','bankssts_right','transversetemporal_right',
                          'superiortemporal_right','temporalpole_right','entorhinal_right',
                          'middletemporal_right','inferiortemporal_right','fusiform_right'],                             
                         ['superiorfrontal_left','frontalpole_left','superiorfrontal_right','frontalpole_right'], 
                         ['caudalmiddlefrontal_left','rostralmiddlefrontal_left','caudalmiddlefrontal_right','rostralmiddlefrontal_right'], 
                         ['parsopercularis_left','parsorbitalis_left','parstriangularis_left','parsopercularis_right','parsorbitalis_right','parstriangularis_right'], 
                         ['medialorbitofrontal_left','medialorbitofrontal_right'], ['lateralorbitofrontal_left','lateralorbitofrontal_right'], 
                         ['precentral_left','paracentral_left','precentral_right','paracentral_right'], ['postcentral_left','postcentral_right'], 
                         ['superiorparietal_left','precuneus_left','superiorparietal_right','precuneus_right','cuneus_left','pericalcarine_left','cuneus_right','pericalcarine_right'], ['inferiorparietal_left','supramarginal_left','inferiorparietal_right','supramarginal_right'], 
                         ['lateraloccipital_left','lateraloccipital_right'],  
                         ['lingual_left','lingual_right'], ['insula_left','insula_right'], ['caudalanteriorcingulate_left','rostralanteriorcingulate_left','caudalanteriorcingulate_right','rostralanteriorcingulate_right'], 
                         ['posteriorcingulate_left','isthmuscingulate_left','posteriorcingulate_right','isthmuscingulate_right'], 
                         ['parahippocampal_left','parahippocampal_right']]

    dic = {}
    
    if hemi == False:
        for idx, i in enumerate(list_cortical):
            dic[i.lower()] = org_cortical_mapping[idx] 
            
    elif hemi == True:
        list_cortical = list_cortical_left + list_cortical_right
        list_mapping = org_cortical_mapping_left + org_cortical_mapping_right
    
        for idx, i in enumerate(list_cortical):
            dic[i.lower()] = list_mapping[idx] 
        
    return dic


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
           
    dic = dict(zip(labels_ordered, T.sequence_model['event_centers'][subtype_labels[subtype]]))
                    
    # flat lost of dict values (single list of DK-labels)
    dk_flat = [x for v in mapped_dict.values() for x in v]
            
    #Match T.biomarker_labels to DK labels
    list_plot = list()
    for key in mapped_dict.keys():
            for item in mapped_dict[key]:
                list_plot.append(dic[key])
                    
    # Dict for dk-label: T.label value
    dic_dk = dict(zip(dk_flat, list_plot))
    
    return dic_dk


def aseg_dict(T, S, subtype_labels = None, subtype = None, hemi = False):
    
    """
    Creates a dictionary, which can be used as input to ggseg.plot_dk() function
    :param T: dataframe from dk_dataframe() function
    :param S: chosen subtype
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
    labels_dict = {num: label for num, label in enumerate(labels)}
    order = T.sequence_model['ordering'][subtype_labels[subtype]]
    
    labels_ordered = []
    for o in order:
        labels_ordered.append(labels_dict[o])    
    
    if hemi == True:
        dic_aseg = dict(zip(labels_ordered, T.sequence_model['event_centers'][subtype_labels[subtype]]))
    else:
        regions= labels_ordered*2

        for l in range(len(T.biomarker_labels)):
            regions[l]='Left-' + regions[l]
        for l in range(len(T.biomarker_labels), 2*len(T.biomarker_labels)):
            regions[l]='Right-' + regions[l]
            
        values = list(T.sequence_model['event_centers'][subtype_labels[subtype]]) + list(T.sequence_model['event_centers'][subtype_labels[subtype]])                                                                   
                        
        dic_aseg = dict(zip(regions, values))
        
        
    return dic_aseg

# ====================================== DK-ATLAS ==================================================================================================================


def dk_regions_2D(T):
    """
    Creates a dictionary of DK-atlas labels grouped into larger regions corresponding to T.biomarker_labels
    :param T: Timeline object from snowphlake
    :return: dictionary, key:value => T.biomarker_labels: [dk labels]
    """   
    org_cortical_mapping_left = [['lh_bankssts_volume','lh_transversetemporal_volume',
                                  'lh_superiortemporal_volume','lh_temporalpole_volume','lh_entorhinal_volume',
                                  'lh_middletemporal_volume','lh_inferiortemporal_volume','lh_fusiform_volume'], 
                                 ['lh_superiorfrontal_volume','lh_frontalpole_volume'], 
                                 ['lh_caudalmiddlefrontal_volume','lh_rostralmiddlefrontal_volume'], 
                                 ['lh_parsopercularis_volume','lh_parsorbitalis_volume','lh_parstriangularis_volume'], 
                                 ['lh_medialorbitofrontal_volume'], ['lh_lateralorbitofrontal_volume'], 
                                 ['lh_precentral_volume','lh_paracentral_volume'], ['lh_postcentral_volume'], 
                                 ['lh_superiorparietal_volume','lh_precuneus_volume'], ['lh_inferiorparietal_volume','lh_supramarginal_volume'], 
                                 ['lh_lateraloccipital_volume'], ['lh_cuneus_volume','lh_pericalcarine_volume'], 
                                 ['lh_lingual_volume'], ['lh_insula_volume'], ['lh_caudalanteriorcingulate_volume','lh_rostralanteriorcingulate_volume'], 
                                 ['lh_posteriorcingulate_volume','lh_isthmuscingulate_volume'], 
                                 ['lh_parahippocampal_volume']]

    list_imaging_cortical_left = ['Temporal_lobe_left','Superior_frontal_gyrus_left',
                                  'Middle_frontal_gyrus_left','Inferior_frontal_gyrus_left', 
                                  'Gyrus_rectus_left','Orbitofrontal_gyri_left','Precentral_gyrus_left',
                                  'Postcentral_gyrus_left','Superior_parietal_gyrus_left', 
                                  'Inferolateral_remainder_of_parietal_lobe_left',
                                  'Lateral_remainder_of_occipital_lobe_left','Cuneus_left','Lingual_gyrus_left', 
                                  'Insula_left','Gyrus_cinguli_anterior_part_left','Gyrus_cinguli_posterior_part_left',
                                  'Parahippocampal_and_ambient_gyri_left']

    org_cortical_mapping_right = [['rh_bankssts_volume','rh_transversetemporal_volume',
                                  'rh_superiortemporal_volume','rh_temporalpole_volume','rh_entorhinal_volume',
                                  'rh_middletemporal_volume','rh_inferiortemporal_volume','rh_fusiform_volume'], 
                                 ['rh_superiorfrontal_volume','rh_frontalpole_volume'], 
                                 ['rh_caudalmiddlefrontal_volume','rh_rostralmiddlefrontal_volume'], 
                                 ['rh_parsopercularis_volume','rh_parsorbitalis_volume','rh_parstriangularis_volume'], 
                                 ['rh_medialorbitofrontal_volume'], ['rh_lateralorbitofrontal_volume'], 
                                 ['rh_precentral_volume','rh_paracentral_volume'], ['rh_postcentral_volume'], 
                                 ['rh_superiorparietal_volume','rh_precuneus_volume'], ['rh_inferiorparietal_volume','rh_supramarginal_volume'], 
                                 ['rh_lateraloccipital_volume'], ['rh_cuneus_volume','rh_pericalcarine_volume'], 
                                 ['rh_lingual_volume'], ['rh_insula_volume'], ['rh_caudalanteriorcingulate_volume','rh_rostralanteriorcingulate_volume'], 
                                 ['rh_posteriorcingulate_volume','rh_isthmuscingulate_volume'], 
                                 ['rh_parahippocampal_volume']]

    list_imaging_cortical_right = ['Temporal_lobe_right',
                                   'Superior_frontal_gyrus_right',
                                  'Middle_frontal_gyrus_right',
                                   'Inferior_frontal_gyrus_right', 
                                  'Gyrus_rectus_right',
                                   'Orbitofrontal_gyri_right',
                                   'Precentral_gyrus_right',
                                  'Postcentral_gyrus_right',
                                   'Superior_parietal_gyrus_right', 
                                  'Inferolateral_remainder_of_parietal_lobe_right',
                                  'Lateral_remainder_of_occipital_lobe_right',
                                   'Cuneus_right',
                                   'Lingual_gyrus_right', 
                                  'Insula_right',
                                   'Gyrus_cinguli_anterior_part_right',
                                   'Gyrus_cinguli_posterior_part_right',
                                  'Parahippocampal_and_ambient_gyri_right']
    
    # DK-labels in left hemisphere grouped into cortical regions corresponding to T.biomarker_labels
    dk_left = [org_cortical_mapping_left[0] + org_cortical_mapping_left[16],
         org_cortical_mapping_left[1] + org_cortical_mapping_left[2]+org_cortical_mapping_left[3]+
         org_cortical_mapping_left[4]+org_cortical_mapping_left[5]+org_cortical_mapping_left[6],
         org_cortical_mapping_left[7]+org_cortical_mapping_left[8]+org_cortical_mapping_left[9],
         org_cortical_mapping_left[10]+org_cortical_mapping_left[11]+org_cortical_mapping_left[12],
         org_cortical_mapping_left[14]+org_cortical_mapping_left[15],
         org_cortical_mapping_left[13]]
    
    # DK-labels in right hemisphere grouped into cortical regions corresponding to T.biomarker_labels
    dk_right = [org_cortical_mapping_right[0] + org_cortical_mapping_right[16],
         org_cortical_mapping_right[1] + org_cortical_mapping_right[2]+org_cortical_mapping_right[3]+
         org_cortical_mapping_right[4]+org_cortical_mapping_right[5]+org_cortical_mapping_right[6],
         org_cortical_mapping_right[7]+org_cortical_mapping_right[8]+org_cortical_mapping_right[9],
         org_cortical_mapping_right[10]+org_cortical_mapping_right[11]+org_cortical_mapping_right[12],
         org_cortical_mapping_right[14]+org_cortical_mapping_right[15],
         org_cortical_mapping_right[13]]
    
    # clean region names
    for l in range(len(dk_left)):
        for i in range(len(dk_left[l])):
            dk_left[l][i]=dk_left[l][i].replace('_volume','_left')
            dk_left[l][i]=dk_left[l][i].replace('lh_','')
    
    for l in range(len(dk_right)):
        for i in range(len(dk_right[l])):
            dk_right[l][i]=dk_right[l][i].replace('_volume','_right')
            dk_right[l][i]=dk_right[l][i].replace('rh_','')
    
    dk = dk_left + dk_right

    
    regions = list(map(lambda x: x.lower(), T.biomarker_labels[12:]))
    
    # final dictionary of key: value pairs corresponding to T.biomarker_label: list(DK-labels)
    dic = dict(zip(regions, dk))
    
    return dic


def dk_dict_2D(T,S, mapped_dict, subtype_labels = None, subtype = None):
    
    """
    Creates a dictionary, which can be used as input to ggseg.plot_dk() function
    :param T: timeline object from snowphlake
    :param S: dictionary from snowphlake
    :param mapped_dict: mapping_dk() output; a dictionary with key: values --> T.biomarker_labels: list(DK-labels)
    :param subtype: name or index of the subtype from subtype_lables (optional)  
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
       
    dic = dict(zip(labels, T.sequence_model['event_centers'][subtype_labels[subtype]]))
                
    # flat lost of dict values (single list of DK-labels)
    dk_flat = [x for v in mapped_dict.values() for x in v]
        
    #Match T.biomarker_labels to DK labels
    list_plot = list()
    for key in mapped_dict.keys():
        for item in mapped_dict[key]:
            list_plot.append(dic[key])
    
    # Dict for dk-label: T.label value
    dic_dk = dict(zip(dk_flat, list_plot))
    
    return dic_dk

# ====================================== ASEG atlas ==================================================================================================================


def aseg_dict_2D(T, S, subtype_labels = None, subtype = None):
    
    """
    Creates a dictionary, which can be used as input to ggseg.plot_dk() function from R ggseg3d package
    :param T: timeline object from snowphlake
    :param S: dictionary from snowphlake
    :param subtype_labels: a list with names of the subtypes (optional)
    :param subtype: name or index of the subtype from subtype_lables (optional)  
    :return: dictionary with scores for each DK region for chosen subtype
    """

    unique_subtypes = np.unique(S['subtypes'][~np.isnan(S['subtypes'])])
    if subtype_labels is None:
        subtype_labels = {f'Subtype {i}': i for i in range(len(unique_subtypes))}
        if subtype is None:
            subtype = next(iter(subtype_labels))
    elif subtype is None:
        subtype = subtype_labels[0]
    
    dic_aseg = dict(zip(T.biomarker_labels, T.sequence_model['event_centers'][subtype_labels[subtype]]))
        
    return dic_aseg




