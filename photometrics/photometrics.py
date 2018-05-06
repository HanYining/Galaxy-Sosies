import pandas as pd
from sklearn.decomposition import PCA
import re
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import numpy as np
import scipy
import matplotlib.pyplot as plt
import seaborn as sns


class photo_processor:

    def __init__(self, filename):
        self.data = pd.read_csv(filename)
        self._mag = [col for col in self.data.columns if re.search("Mag", col)]
        self._refer = [col for col in self._mag if re.search("Mag_g", col)]
        self._mag = [col for col in self._mag if col not in self._refer]
        self.redshift = self.data["specz"]

    def _take_ratio(self):

        self.cnt = 0
        self.refcnt = -1

        for col in self._mag:
            self.cnt += 1
            if self.cnt % 4 == 1:
                self.refcnt += 1

            self.data[col] = self.data[col]/self.data[self._refer[self.refcnt]]

        self.data.index = self.data["specObjID"]
        self.data = self.data[self._mag]

        pass

    def bindt(self, num_bins=10):
        """bin the original dataset to produce binned dataframe on redshift"""
        self.bin_data = []
        quantiles = [i*(float(1)/num_bins) for i in range(num_bins+1)]
        quantiles_redshift = self.redshift.quantile(
            quantiles, interpolation="linear")

        for i in range(len(quantiles_redshift)-1):
            ind = (quantiles_redshift.iloc[i] <= self.redshift) \
                & (self.redshift < quantiles_redshift.iloc[i+1])
            ind = ind.reshape(-1, )
            self.bin_data.append(self.data.loc[ind])

        pass

    def photometrics_pca(self, threshold=0.75):
        """
        Implement a principal components analysis on the
        emission fluxes, preserving threshold % amount of
        variance in the original dataset.
        remaining principal directions are labeled as
        "flux1", "flux2", etc
        """
        for i, bin in enumerate(self.bin_data):

            # pca on flux
            ind = bin.index
            scaler = StandardScaler().fit(bin)
            data = scaler.transform(bin)

            j = 1
            pca = PCA(n_components=j)
            pca.fit(data)
            while sum(pca.explained_variance_ratio_) < threshold:
                j += 1
                pca = PCA(n_components=j)
                pca.fit(data)

            flux_pca = pca.fit_transform(data)
            photometrics_column = ["photo" + str(k+1) for k in range(j)]
            flux_df = pd.DataFrame(data=flux_pca,
                                   columns=photometrics_column)
            flux_df.index = ind
            self.bin_data[i] = flux_df

        pass
