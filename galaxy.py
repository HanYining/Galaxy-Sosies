import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
import numpy as np
import re


class galaxy_processor:

    def __init__(self, filename):
        self._data = pd.read_csv(filename)
        self._data = self._data[[column for column in self._data.columns if
                                not re.search("Mag", column)]]

        self._emissions = ['Flux_HeII_3203', 'Flux_NeV_3345', 'Flux_NeV_3425',
                           'Flux_OII_3726', 'Flux_OII_3728', 'Flux_NeIII_3868',
                           'Flux_NeIII_3967', 'Flux_H5_3889', 'Flux_He_3970',
                           'Flux_Hd_4101', 'Flux_Hg_4340', 'Flux_OIII_4363',
                           'Flux_HeII_4685', 'Flux_ArIV_4711',
                           'Flux_ArIV_4740', 'Flux_Hb_4861',
                           'Flux_OIII_4958', 'Flux_OIII_5006',
                           'Flux_NI_5197', 'Flux_NI_5200', 'Flux_HeI_5875',
                           'Flux_OI_6300', 'Flux_OI_6363', 'Flux_NII_6547',
                           'Flux_Ha_6562', 'Flux_NII_6583', 'Flux_SII_6716',
                           'Flux_SII_6730']

        self._photometrics = [col for col in self._data.columns
                              if col not in self._emissions and
                              col not in ["specObjID", "objid", "ra", "dec",
                                          "photoz", "specz", "type"]]

        self._identifications = ["specObjID", "objid"]

        self._objective = ["specz"]

        self._bin_data = []

    def zero_est(self, data):
        """
        a brief function to check the non-emission estimates of emission lines
        output a percentage of non-emission estimates
        """
        invalid_map = {}
        for column in self._emissions:
            invalid_map.update(
                {column: round(100*np.sum(data[column] == 0)
                               / data.shape[0], 2)})
        return invalid_map

    def negative_est(self, data):
        """a function to check the error estimates"""

        missing_map = {}
        print("\n")
        for column in self._emissions:
            missing_map.update(
                {column: round(100*np.sum(data[column] < 0)
                               / data.shape[0], 2)})
        return missing_map

    def bin_data(self, num_bins=10):
        """bin the original dataset to produce binned dataframe on redshift"""

        quantiles = [i*(float(1)/num_bins) for i in range(num_bins+1)]
        quantiles_redshift = self._data[self._objective].quantile(
            quantiles, interpolation="linear")

        for i in range(len(quantiles_redshift)-1):
            self._bin_data.append(self._data.loc[
                (quantiles_redshift.iloc[i][0] <= self._data["specz"]) &
                (self._data["specz"] < quantiles_redshift.iloc[i+1][0])
            ])
        pass

    def filter_valid(self):
        """filter through the observations that has negative"""
        for i, bin in enumerate(self._bin_data):
            filter = np.all([bin[col] >= 0 for col in self._emissions], axis=0)
            self._bin_data[i] = bin[filter]
        pass

    def filter_under(self, missing_percent=50):
        """filter through every features which has less than 30% of 0s"""
        for i, bin in enumerate(self._bin_data):
            zero_estimation = self.zero_est(bin)
            ok_feature = [col for col, perc in zero_estimation.items()
                          if perc <= missing_percent]
            self._bin_data[i] = bin[self._identifications +
                                    self._photometrics + ok_feature]
            print("bin{0}: {1}".format(i, ok_feature))

        pass

    def normalize_flux(self):
        """
        Normalizing the emission flux by the emission
        Flux which has the highest average.
        Recommended but does not really make sense to me
        """

        for i, bin in enumerate(self._bin_data):
            emissions = [string for string in bin.columns
                         if re.search("Flux", string)]
            means = np.mean(bin[emissions], axis=0)
            column = means.idxmax()
            normalizing_column = bin[column]
            normalizing_column[normalizing_column == 0] = \
                np.mean(normalizing_column)
            print("The reference Flux in bin:{0}: {1}".
                  format(i+1, column))
            new_data = bin[[string for string in emissions
                            if string != column]]
            origin_data = bin[[string for string in bin.columns
                               if string not in new_data.columns and
                               string != column]]
            normalizing_column = normalizing_column.reshape(-1, 1)
            new_data = np.divide(new_data, normalizing_column)
            self._bin_data[i] = pd.concat([new_data, origin_data],
                                          axis=1)

        pass

    def flux_pca(self, threshold=0.85):
        """
        Implement a principal components analysis on the
        emission fluxes, preserving threshold % amount of
        variance in the original dataset.
        remaining principal directions are labeled as
        "flux1", "flux2", etc
        """
        for i, bin in enumerate(self._bin_data):

            # pca on flux
            emissions = [string for string in bin.columns
                         if re.search("Flux", string)]
            pca = PCA()
            pca.fit(bin[emissions])
            j = 1
            while sum(pca.explained_variance_ratio_[:j]) < threshold:
                j += 1
            pca = PCA(n_components=j)
            flux_pca = pca.fit_transform(bin[emissions])
            flux_columns = ["flux" + str(k+1) for k in range(j)]
            flux_df = pd.DataFrame(data=flux_pca, columns=flux_columns)
            origin_data = bin[[col for col in bin.columns
                               if col not in emissions]]
            flux_df.index = origin_data.index
            self._bin_data[i] = pd.concat([origin_data, flux_df], axis=1)

        pass

    def generate_cluster(self):
        """
        not sure if we want to normalize or standardize on the
        photometrics side though.

        DBSCAN implementation here.
        This function generate the cluster in each bin
        based on the photometrics feature.

        reference:
        http://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
        http://www.aaai.org/Papers/KDD/1996/KDD96-037.pdf
        """

        for i, bin in enumerate(self._bin_data):

            dbscan_cluster = DBSCAN()
            labels = dbscan_cluster.fit_predict(bin[self._photometrics])
            labels = labels.reshape(-1, 1)
            label = pd.DataFrame(data=labels, columns=["group"])
            label.index = bin.index
            new_data = pd.concat([bin, label], axis=1)
            new_data = new_data[new_data["group"] != -1]
            print("distinct groups in bin{0} :{1}".
                  format(i+1, len(np.unique(label))-1))

            self._bin_data[i] = new_data

        pass




glx = galaxy_processor("summary_yinhan.csv")
glx.bin_data()
glx.filter_valid()
glx.filter_under()
glx.normalize_flux()
glx.flux_pca()
glx.generate_cluster()

# now the following is a list of each bins
# every bin contains a dataframe which represents the dimension reduced
# galaxies.
# also the group column represents the group effects from photometrics.

