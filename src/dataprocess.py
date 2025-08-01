import os
import sys
import numpy as np
import pandas as pd
new_path = os.path.abspath(os.path.join(__file__, "../../"))
sys.path[1] = new_path
feature_path = sys.path[1] + '/Result/figures/Feature/'
data_path = sys.path[1] + '/Result/data/'
Norm_path=sys.path[1] + '/Result/data/'
def split_data(analysis_level, groups,data_type):
    data = Norm_path + f"{data_type}/Raw/{analysis_level}/Raw_feature_finally.csv"
    feature = pd.read_csv(data, sep=',', index_col=0)
    meta = pd.read_csv(data_path+f"{data_type}/Raw_meta_finally.csv")
    merged_df = pd.merge(feature, meta[["Group", "SampleID", "Study"]], left_index=True, right_on='SampleID', how="inner")
    merged_df.set_index('SampleID', inplace=True)
    return merged_df
def select_feature(analysis_level, group_name, feature_name,data_type):
    feature_dir = (feature_path+f"{data_type}/{group_name}/{analysis_level}/All_features_tools/All_features_adjP.csv")
    feature = pd.read_csv(feature_dir, index_col=0)
    feature_select = feature[feature[feature_name]<0.05].index.tolist()
    return feature_select
def balance_samples(meta_feature, group1, group2):
    grouped = meta_feature.groupby("Study")
    selected_samples = pd.DataFrame()
    for name, group in grouped:
        ada_samples = group[group["Group"] == group1]
        ctr_samples = group[group["Group"] == group2]
        min_count = min(len(ada_samples), len(ctr_samples))
        if min_count > 0:
            ada_samples = ada_samples.sample(n=min_count, replace=False)
            ctr_samples = ctr_samples.sample(n=min_count, replace=False)
            selected_samples = selected_samples.append(ada_samples)
            selected_samples = selected_samples.append(ctr_samples)
    selected_samples.reset_index(drop=False, inplace=True)
    selected_samples.set_index('Sample_ID', inplace=True)
    return selected_samples

def change_group(meta_feature,group_name):
    """
      Convert group labels to binary.
      :param meta_feature: Metadata dataframe with group column.
      :param group_name: Name of the target group.
      :return: Updated metadata dataframe with binary group labels.
      """
    if group_name == 'CTR_ADA':
        meta_feature['Group'] = meta_feature['Group'].apply(lambda x: 1 if x == "ADA" else 0)
    elif group_name == "CTR_CRC":
        meta_feature['Group'] = meta_feature['Group'].apply(lambda x: 1 if x == "CRC" else 0)
    else:
        meta_feature['Group'] = meta_feature['Group'].apply(lambda x: 1 if x == "CRC" else 0)
    return meta_feature

