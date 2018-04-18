import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn import preprocessing
import seaborn as sns
import numpy as np
import random


def invalid_est(df):

    """
    a brief function to check the invalid estimates of emission lines
    output a percentage of missing estimates
    """

    print("\n")
    for column in df.columns:
        print("{0}:{1}%".format(
            column, round(100*np.sum(df[column] == 0)/df.shape[0], 2)))


def pie_chart(df):

    """
    This function returns a random galaxy pie chart plot from all the
    galaxy data frames.

    This takes cares of the case
    (etc that there are multiple NeIII emission lines)
    in the data frame and add them up to give a total estimate of the
    particular component.
    """

    idx = random.randint(0, df.shape[0])
    bool_index = df.loc[idx] != 0
    label = [column.split("_")[1] for column, ind in
             zip(df.columns, bool_index) if ind]
    values = df.loc[idx][bool_index]

    distinct_label = set(label)
    distinct_value = []
    for lab in distinct_label:
        current_value = 0
        for i, lab2 in enumerate(label):
            if lab2 == lab:
                current_value += values[i]
        distinct_value.append(current_value)

    plt.pie(distinct_value, labels=distinct_label,
            autopct='%1.1f%%', shadow=True, startangle=140)

    plt.show()


data = pd.read_csv("summary_yinhan.csv")
plt.hist(data["Flux_NI_5200"])
print(sum(data["Flux_He_3970"] <= 0)/data.shape[0])
plt.show()

data.replace(to_replace=-9999, value=0, inplace=True)
emission_lines = ['Flux_HeII_3203', 'Flux_NeV_3345', 'Flux_NeV_3425',
                  'Flux_OII_3726', 'Flux_OII_3728', 'Flux_NeIII_3868',
                  'Flux_NeIII_3967', 'Flux_H5_3889', 'Flux_He_3970',
                  'Flux_Hd_4101', 'Flux_Hg_4340', 'Flux_OIII_4363',
                  'Flux_HeII_4685', 'Flux_ArIV_4711', 'Flux_ArIV_4740',
                  'Flux_Hb_4861', 'Flux_OIII_4958', 'Flux_OIII_5006',
                  'Flux_NI_5197', 'Flux_NI_5200', 'Flux_HeI_5875',
                  'Flux_OI_6300', 'Flux_OI_6363', 'Flux_NII_6547',
                  'Flux_Ha_6562', 'Flux_NII_6583', 'Flux_SII_6716',
                  'Flux_SII_6730']

emission_data = data[emission_lines]

# PCA on the original whole dataset
# notice some of those emission lines has extreme flux
pie_chart(emission_data)
invalid_est(emission_data)
scaler = preprocessing.StandardScaler().fit(emission_data)
emission_data_std = scaler.transform(emission_data)
pca = PCA(n_components=6)
pca.fit(emission_data_std)
print(sum(pca.explained_variance_ratio_))

sns.distplot(emission_data["Flux_NeIII_3967"])

# filter through the galaxies using 95% quantile on ALL features
# this left us with only about a half of the original galaxies.
# not sure if this makes sense to you.

filter_quantile = np.all([emission_data[col] < emission_data[col].quantile(.95)
                          for col in emission_lines], axis=0)
filtered_galaxy = emission_data[filter_quantile]
invalid_est(filtered_galaxy)
scaler_filtered = preprocessing.StandardScaler().fit(filtered_galaxy)
filtered_galaxy_std = scaler_filtered.transform(filtered_galaxy)
pca = PCA(n_components=6)
pca.fit(filtered_galaxy_std)
print(sum(pca.explained_variance_ratio_))

sns.distplot(filtered_galaxy["Flux_HeII_4685"])
plt.show()
