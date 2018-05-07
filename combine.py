from ratio import ratio
from photometrics import photometrics
from flux import flux
import pandas as pd
import numpy as np
import scipy
import re
from sklearn.cluster import DBSCAN


# combine the ratio of AB
rat = ratio.galaxy_ratio("summary_yinhan.csv")
rat.average()
rat.bin_ratio(ratio_bins=10)

# combine the photometrics features
photo = photometrics.photo_processor("summary_yinhan.csv")
photo._take_ratio()
photo.bindt(10)
photo.photometrics_pca()

# combine the fluxes reduced features
flx = flux.flux_operation("summary_yinhan.csv")
flx.bin_data(10)
flx.filter_valid()
flx.filter_anomaly()
flx.normalize_flux()
flx.flux_pca()


def combine_redshift_bins(lst):
    for i, data in enumerate(lst):
        data["redshift_group"] = pd.Series(np.repeat(i+1, data.shape[0]),
                                           data.index)
        lst[i] = data
    return lst


flx_data = pd.concat(combine_redshift_bins(flx.bin_data), axis=0)
photo_data = pd.concat(combine_redshift_bins(photo.bin_data), axis=0)
photo_data = photo_data.drop(["redshift_group"], axis=1)
ratio_data = rat.data
data = pd.concat([photo_data, ratio_data, flx_data], join="inner", axis=1)


# partition dataset

def _get_distance_dist(data, threshold):
    """
    get the pairwise distance from the photometrics features
    here because the size of the data, pairwise operations can
    quickly blow up, i choose to randomly sample a subset and estimate
    the distance distribution.

    this function is called by the cluster generator thus we could have
    a sense on how close average points are and we can estimate a good
    value for the epsilon in dbscan.
    """

    epsilon_distance = []
    dist = scipy.spatial.distance.pdist(data, metric='euclidean')
    epsilon_distance.append(np.percentile(dist, [threshold]))

    return epsilon_distance


def generate_cluster(data, threshold):
    """
    not sure if we want to normalize or standardize on the
    photometrics side though.

    dbscan implementation here.
    this function generate the cluster in each bin
    based on the photometrics feature.

    reference:
    http://scikit-learn.org/stable/modules/generated/sklearn.cluster.dbscan.html
    http://www.aaai.org/papers/kdd/1996/kdd96-037.pdf
    """

    epsilon_distance = _get_distance_dist(data, threshold)[0][0]

    dbscan_cluster = DBSCAN(eps=epsilon_distance)
    labels = dbscan_cluster.fit_predict(data)
    labels = labels.reshape(-1, 1)
    label = pd.DataFrame(data=labels, columns=["photo_group"])
    label.index = data.index
    print("distinct groups {0}".format(len(np.unique(label))-1))

    return label.loc[np.all([label["photo_group"] != -1,
                             label["photo_group"] != 0], axis=0)]


def data_photo_group(data):

    ratio_group_num = int(max(data["ratio_group"]))
    redshift_group_num = int(max(data["redshift_group"]))

    label_list = []

    counter = 0
    for rat in range(1, ratio_group_num+1):
        for red in range(1, redshift_group_num+1):
            ind = np.all([data["ratio_group"] == rat,
                          data["redshift_group"] == red], axis=0)
            new_data = data.loc[ind]
            if new_data.shape[0] <= 1:
                continue
            photo_data = new_data[["photo1", "photo2"]]
            label = generate_cluster(photo_data, 1)
            label = label+counter
            counter += label["photo_group"].unique().shape[0]
            label_list.append(label)

    return pd.concat([data, pd.concat(label_list, axis=0)],
                     axis=1, join="inner")


photo_group = data_photo_group(data)


def filter_groups(data):

    crude_group_num = int(max(data["photo_group"]))

    lst = []

    for group in range(1, crude_group_num+1):

        ind = data["photo_group"] == group
        columns = [col for col in data.columns if re.search("flux", col)]
        new_data = data.loc[ind][columns]
        new_data = new_data.dropna(axis=1, how="all")

        mean_group = np.mean(new_data, 0)
        dist = np.mean((new_data - mean_group)**2, axis=1)

        threshold_dist = np.percentile(dist, 60)
        sosie = new_data[np.mean((new_data - mean_group)**2, axis=1)
                         < threshold_dist]

        if sosie.shape[0] <= 1:
            continue

        lst.append(sosie)

    data = pd.concat([pd.concat(lst), data[["photo_group", "ratio_group",
                                           "redshift_group"]]],
                     axis=1, join="inner")
    data = data.rename({"photo_group":"sosie_group"}, axis="columns")
    columns = [col for col in data.columns if not re.search("flux", col)]

    return data[columns]


refined_group = filter_groups(photo_group)
data = pd.read_csv("summary_yinhan.csv")
data.index = data["objid"]
final_data = pd.concat([refined_group, data[["cModelMag_g", "specz"]]],
                       axis=1, join="inner")

final_data["lower_redshift"] = final_data["specz"]
final_data["higher_redshift"] = final_data["specz"]

num_bins = 10
quantiles_redshift = list(np.linspace(0, 0.9, 10)) + [100]

for i in range(num_bins):
    final_data.loc[final_data["redshift_group"] == i+1,
                   "lower_redshift"] = quantiles_redshift[i]
    final_data.loc[final_data["redshift_group"] == i+1,
                   "higher_redshift"] = quantiles_redshift[i+1]

final_data["ratio_group"] = final_data.ratio_group.astype(int)

final_data.to_csv("sosie_list.txt", sep="\t")


def _calculate_distance(df):
    """
    this function takes in a pandas dataframe with columns representing the
    luminosity(by petromag_g) and their redshift specz
    for which we are thinking about give a specific minimum ratio

    here because of the relative scale of the petromag_g
    there is not so much difference in terms of the absolute value
    and there is a nasty reference line i have to choose.
    here i set it to 14
    """

    reference_magnitude = 0

    magnitude = []
    specz = []
    df = df.sort_values(by=["cModelMag_g", "specz"])
    magnitude.append((df.iloc[-1]["cModelMag_g"] - reference_magnitude)
                     / (df.iloc[0]["cModelMag_g"] - reference_magnitude))
    specz.append(df.iloc[-1]["specz"]/df.iloc[0]["specz"])

    return [sum(magnitude)/len(magnitude), sum(specz)/len(specz), df.shape[0]]


def point_generater(sosies_data):

    res = []
    for group, df in sosies_data.groupby('sosie_group'):
        res.append(_calculate_distance(df))

    data = pd.DataFrame(data=res, columns=["distance", "redshift", "weight"])
    data["distance"] = np.sqrt(data["distance"])

    return data
