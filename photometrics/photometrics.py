import pandas as pd
from sklearn.decomposition import PCA
import re
import numpy as np
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt


class photo_processor:

    def __init__(self, filename):
        self.data = pd.read_csv(filename)
        self._mag = [col for col in self.data.columns if re.search("Mag", col)]
        self._refer = [col for col in self._mag if re.search("Mag_g", col)]
        self._mag = [col for col in self._mag if col not in self._refer]
        self.redshift = self.data["specz"]

        self.data.index = self.data["objid"]
        self.data = self.data[self._mag + self._refer]

    def _take_ratio(self):

        self.cnt = 0
        self.refcnt = -1

        for col in self._mag:
            self.cnt += 1
            if self.cnt % 4 == 1:
                self.refcnt += 1

            self.data[col] = self.data[col]/self.data[self._refer[self.refcnt]]
        self.data = self.data[self._mag]
        pass

    def bindt(self, num_bins=10):
        """bin the original dataset to produce binned dataframe on redshift"""
        self.bin_data = []
        quantiles_redshift = list(np.linspace(0, 0.9, 10)) + [100]

        for i in range(len(quantiles_redshift)-1):
            ind = (quantiles_redshift[i] <= self.redshift) \
                & (self.redshift < quantiles_redshift[i+1])
            ind = ind.values.reshape(-1, )
            self.bin_data.append(self.data.loc[ind])

        pass

    def photometrics_pca(self):
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
            pca = PCA(n_components=2)
            pca.fit(data)

            flux_pca = pca.fit_transform(data)
            photometrics_column = ["photo" + str(k+1) for k in range(2)]
            flux_df = pd.DataFrame(data=flux_pca,
                                   columns=photometrics_column)
            flux_df.index = ind
            self.bin_data[i] = flux_df

        pass


if __name__ == "__main__":
    # some basic edas
    ph = photo_processor("summary_yinhan.csv")
    data = ph.data
    sub_sample = data.sample(10000)

    # it is clear from those pictures, that the shifts among the features
    # themselves is dominated by the distance.
    # notice here the wavelength rank of the corresponding mode for the
    # distribution is u/g/r/i/z

    sns.distplot(sub_sample["cModelMag_u"], hist=False, label="u")
    sns.distplot(sub_sample["cModelMag_g"], hist=False, label="g")
    sns.distplot(sub_sample["cModelMag_r"], hist=False, label="r")
    sns.distplot(sub_sample["cModelMag_i"], hist=False, label="i")
    sns.distplot(sub_sample["cModelMag_z"], hist=False, label="z")
    plt.legend(loc='upper right')
    plt.xlabel("Magnitude")
    plt.title("Density Estimation Plot by Filter")
    plt.show()

    # this picture shows that there are a lot of correlation between
    # different magnitudes.

    sns.pairplot(sub_sample)
    plt.show()

    # why is this?
    # the first thought should be that because of the distance for each
    # galaxy is different. Far away galaxies naturally have bigger magnitude
    # appear to be dimmer.
    # first decision: take the ratio between the measurements.

    ph._take_ratio()
    data2 = ph.data
    sub_sample2 = data2.sample(10000)
    sns.pairplot(sub_sample2)
    plt.show()

    # now even if the features are invariant to the change of distance.
    # the situation gets better but we still see a very strong correlation

    # now we turn to some statistical methods to lower the dimension
    # and get rid of the correlation between features

    ph.bindt()
    ph.photometrics_pca()

    # simply use a random bin
    data3 = ph.bin_data[4]
    sub_sample3 = data3.sample(10000)
    sns.pairplot(sub_sample3)
    plt.show()

    # now it seems rather safe, the features are no longer highly
    # correlated.
