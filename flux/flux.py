import pandas as pd
import numpy as np
import re
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class flux_operation:

    def __init__(self, filename):
        self.data = pd.read_csv(filename)
        self.redshift = self.data["specz"]
        columns = [col for col in self.data.columns if re.search("Flux", col)]
        self.data.index = self.data["objid"]
        self.data = self.data[columns]

    def zero_est(self, data):
        """
        a brief function to check the non-emission estimates of emission lines
        output a percentage of non-emission estimates
        """
        invalid_map = {}
        for column in data.columns:
            invalid_map.update(
                {column: round(100*np.sum(data[column] == 0)
                               / data.shape[0], 2)})
        return invalid_map

    def negative_est(self, data):
        """a function to check the error estimates"""

        missing_map = {}
        print("\n")
        for column in data.columns:
            missing_map.update(
                {column: round(100*np.sum(data[column] < 0)
                               / data.shape[0], 2)})
        return missing_map

    def bin_data(self, num_bins=10):
        """bin the original dataset to produce binned dataframe on redshift"""

        self.bin_data = []
        quantiles_redshift = list(np.linspace(0, 0.9, 10)) + [100]
        for i in range(len(quantiles_redshift)-1):
            ind = (quantiles_redshift[i] <= self.redshift) \
                & (self.redshift < quantiles_redshift[i+1])
            ind = ind.values.reshape(-1, )
            self.bin_data.append(self.data.loc[ind])

        pass

    def filter_valid(self):
        """filter through the observations that has negative"""
        for i, bin in enumerate(self.bin_data):
            filter = np.all([bin[col] >= 0 for col in bin.columns], axis=0)
            self.bin_data[i] = bin[filter]
        pass

    def filter_anomaly(self):
        for i, bin in enumerate(self.bin_data):
            ind = np.all([np.any([np.percentile(bin[col], 98) > bin[col],
                                  bin[col] == 0], axis=0)
                          for col in bin.columns], axis=0)
            self.bin_data[i] = bin[ind]
        pass

    def normalize_flux(self):
        """
        normalizing the emission flux by the emission
        flux which has the highest average.
        recommended but does not really make sense to me
        """

        for i, bin in enumerate(self.bin_data):
            means = np.mean(bin == 0, axis=0)
            column = means.idxmin()
            bin = bin.loc[bin[column] != 0]
            bin = bin.loc[np.percentile(bin[column], 1) < bin[column]]
            index = bin.index

            normalizing_column = bin[column]

            print("the reference flux in bin:{0}: {1}".
                  format(i+1, column))
            columns = [string for string in bin.columns if string != column]
            new_data = bin[columns]
            normalizing_column = normalizing_column.reshape(-1, 1)
            new_data = np.divide(new_data, normalizing_column)
            self.bin_data[i] = pd.DataFrame(new_data, index=index,
                                            columns=columns)

        pass

    def flux_pca(self, threshold=0.75):
        """
        implement a principal components analysis on the
        emission fluxes, preserving threshold % amount of
        variance in the original dataset.
        remaining principal directions are labeled as
        "flux1", "flux2", etc
        """
        for i, bin in enumerate(self.bin_data):

            # pca on flux
            data = bin
            scaler = StandardScaler().fit(data)
            data = scaler.transform(data)

            j = 1
            pca = PCA(n_components=j)
            pca.fit(data)
            while sum(pca.explained_variance_ratio_) < threshold:
                j += 1
                pca = PCA(n_components=j)
                pca.fit(data)

            flux_pca = pca.fit_transform(data)
            flux_columns = ["flux" + str(k+1) for k in range(j)]
            flux_df = pd.DataFrame(data=flux_pca, index=bin.index,
                                   columns=flux_columns)
            self.bin_data[i] = flux_df

        pass


if __name__ == "__main__":

    data = pd.read_csv("summary_yinhan.csv")
    columns = [col for col in data.columns if re.search("Flux", col)]
    data['redshift_group'] = np.nan
    data["redshift_group"] = data["specz"]
    data.index = data["objid"]
    data = data[columns + ["specz", "redshift_group"]]

    quantiles_redshift = list(np.linspace(0, 0.9, 10)) + [100]
    for i in range(len(quantiles_redshift)-1):
        ind = (quantiles_redshift[i] <= data["specz"]) \
            & (data["specz"] <= quantiles_redshift[i+1])
        ind = ind.values.reshape(-1, )
        data.loc[ind, "redshift_group"] = i+1

    # make a table for the percentage of zeros:

    res = []
    for group, sub_data in data.groupby("redshift_group"):
        res.append(pd.DataFrame(np.mean(sub_data[columns] != 0, axis=0),
                                columns=["Bin"+str(int(group))]))

    res_df = pd.concat(res, axis=1, join="inner")
