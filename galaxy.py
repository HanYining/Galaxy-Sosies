import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn import preprocessing
import seaborn as sns
import numpy as np
import re


class galaxy_processor:

    def __init__(self, filename):
        self._data = pd.read_csv(filename)

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
            print("{0}:{1}%".format(
                column, round(100*np.sum(data[column] == 0)
                              / data.shape[0], 2)))
            invalid_map.update(
                {column: round(100*np.sum(data[column] == 0)
                               / data.shape[0], 2)})
        return invalid_map

    def negative_est(self, data):
        """a function to check the error estimates"""

        missing_map = {}
        print("\n")
        for column in self._emissions:
            print("{0}:{1}%".format(
                column, round(100*np.sum(data[column] < 0)
                              / data.shape[0], 2)))
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

    def post_normalization_pca(self, threshold=0.85):
        """
        This function do a separate pca on flux emission lines and
        photometrics data.
        They are separately implemented and postnormalization is
        applied.
        Also the number of PC vectors chosen on the two sides is by
        using the least number of PC vectors that achieved at least
        the threshold amount of variance explained.
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
            pca = PCA(n_components=j+1)
            flux_pca = pca.fit_transform(bin[emissions])
            flux_pca = flux_pca[:, 1:]/flux_pca[:, 0].reshape(-1, 1)
            flux_columns = ["flux" + str(k) for k in range(j)]
            flux_df = pd.DataFrame(data=flux_pca, columns=flux_columns)

            # pca on photometrics
            pca = PCA()
            pca.fit(bin[self._photometrics])
            j = 1
            while sum(pca.explained_variance_ratio_[:j]) < threshold:
                j += 1
            pca = PCA(n_components=j+1)
            eval_photometrics = pca.fit_transform(bin[self._photometrics])
            eval_photometrics = eval_photometrics[:, 1:] \
                / eval_photometrics[:, 0].reshape(-1, 1)
            photometrics_columns = ["photo" + str(k) for k in range(j)]
            photo_df = pd.DataFrame(data=eval_photometrics,
                                    columns=photometrics_columns)

            iden = bin[self._identifications]

            self._bin_data[i] = pd.concat([iden, flux_df, photo_df])
        pass


glx = galaxy_processor("summary_yinhan.csv")
glx.bin_data()
glx.filter_valid()
glx.filter_under()
glx.post_normalization_pca()

# now the following is a list of each bins
# every bin contains a dataframe which represents the dimension reduced
# galaxies.
# only contains about 1-3 PC from flux and 1-4 features from photometrics
# remaining features look like this on each bin
# ["specobjid", "flux0", "flux1", "photo1", "photo2"]
# notice the features share the same name across bins do not share the
# same meaning.
glx._bin_data
