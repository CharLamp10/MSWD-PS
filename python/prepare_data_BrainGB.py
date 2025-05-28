from scipy.io import loadmat
import numpy as np
import os
import pandas as pd


casee = 'MSWD_'
path = r'C:\Users\100063082\Desktop\MSWD_PS_revision'

rois_path = os.path.join(path, casee + 'rois_ho')
corrs_path = os.path.join(path, casee + 'rois_ho_corr')
rois_dir = os.listdir(rois_path)
corrs_dir = os.listdir(corrs_path)

temp = []
for roi_dir in rois_dir:
    if '.mat' in roi_dir:
        temp.append(roi_dir)
rois_dir = temp

phenotypic = pd.read_csv(r'C:\Users\100063082\Desktop\prepare_ADHD200_ABIDE_phenotypics\Phenotypic_V1_0b_preprocessed1.csv')
ids = phenotypic["SUB_ID"].values
groups = phenotypic["DX_GROUP"].values
age_at_scan = phenotypic['AGE_AT_SCAN'].values
sex = phenotypic['SEX'].values-1

timeseries = np.zeros((1074,111,116))
Label = np.zeros((1074))
corr = np.zeros((1074,111,111))
pcorr = np.zeros((1074,111,111))
ages = np.zeros((1074))
gender = np.zeros((1074))
site = []
for i,(file_roi, file_corr) in enumerate(zip(rois_dir,corrs_dir)):
    pos = file_roi.find('00')
    sub_id = file_roi[pos+2:pos+7]
    pos_sub = np.where(ids == int(sub_id))
    label = groups[pos_sub]
    if label == 2:
        label = 0
    Label[i] = label
    site.append(file_roi[:pos-1])
    if 'MSWD' not in casee:
        time = loadmat(os.path.join(rois_path,file_roi))["data"]
        cor = loadmat(os.path.join(corrs_path,file_corr))["corr"]
    else:
        time = loadmat(os.path.join(rois_path,file_roi))["rec"]
        cor = loadmat(os.path.join(corrs_path,file_corr))["corr_dec"]
    timeseries[i] = time[:116,:].T
    corr[i] = cor
    pcorr[i] = cor
    ages[i] = age_at_scan[pos_sub]
    gender[i] = sex[pos_sub]

df = pd.DataFrame({"age": ages, "gender": gender, "Label": Label})    
df.to_csv("age_gender.csv", index=False)

site = np.array(site)
abide = {}
abide['timeseries'] = timeseries
abide['label'] = Label
abide['corr'] = corr
abide['pcorr'] = pcorr
abide['site'] = site
np.save('abide_MSWD.npy',abide)
