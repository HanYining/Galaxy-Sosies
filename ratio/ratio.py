import pandas as pd
import numpy as np
import re


class galaxy_ratio:

    def __init__(self, filename):

        data = pd.read_csv(filename)
        data.index = data["specObjID"]
        columns = [col for col in data.columns if re.search("AB_", col)]

        self.redshift = data["specz"]
        self.data = data[columns]

    def average(self):

        index = self.data.index
        new_data = np.mean(self.data, axis=1)
        self.data = pd.DataFrame(new_data, index=index, columns=["abratio"])

        pass

    def bin_ratio(self, ratio_bins=10):
        """bin the original dataset to produce binned dataframe on redshift"""
        quantiles = [i*(float(1)/ratio_bins) for i in range(ratio_bins+1)]
        quantiles_abratio = self.data["abratio"].quantile(
            quantiles, interpolation="linear")

        self.data['ratio_group'] = pd.Series(
            np.random.randn(self.data.shape[0]), index=self.data.index)

        for i in range(len(quantiles_abratio)-1):
            ind = (quantiles_abratio.iloc[i] <= self.data["abratio"]) \
                & (self.data["abratio"] < quantiles_abratio.iloc[i+1])
            ind = ind.reshape(-1, )
            self.data["ratio_group"].loc[ind] = np.repeat(i+1, sum(ind))
        self.data = self.data.drop(["abratio"], axis=1)

        pass
