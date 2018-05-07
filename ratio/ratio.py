import pandas as pd
import numpy as np
import re
import seaborn as sns
import matplotlib.pyplot as plt


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


if __name__ == "__main__":

    # examine the relationship between the redshift and ab ratio
    # does the ratio's distribution look similar across redshift bins.

    data = pd.read_csv("summary_yinhan.csv")
    data.index = data["specObjID"]
    columns = [col for col in data.columns if re.search("AB_", col)]
    redshift = data["specz"]
    data2 = pd.DataFrame(np.mean(data[columns], axis=1), columns=["ABratio"])
    data2.index = data.index
    data2 = pd.concat([data2, redshift], axis=1)
    data2["redshift_group"] = data2["ABratio"]

    num_bins = 10
    quantiles = [i*(float(1)/num_bins) for i in range(num_bins+1)]
    quantiles_redshift = data2["specz"].quantile(
        quantiles, interpolation="linear")

    for i in range(len(quantiles_redshift)-1):
        ind = (quantiles_redshift.iloc[i] <= data2["specz"]) \
            & (data2["specz"] < quantiles_redshift.iloc[i+1])
        ind = ind.values.reshape(-1, )
        data2.loc[ind, "redshift_group"] = i

    data = data2.sample(20000)

    # the take away message here is that we should
    # only compare a/b ratios inside each redshift bin
    # otherwise the a/b ratio is skewed by redshift
    # This justify our choice of restraining on a/b ratio
    # after redshift

    sns.kdeplot(data.loc[data["redshift_group"] == 1]["ABratio"],
                shade=True, label="bin1")
    sns.kdeplot(data.loc[data["redshift_group"] == 3]["ABratio"],
                shade=True, label="bin3")
    sns.kdeplot(data.loc[data["redshift_group"] == 5]["ABratio"],
                shade=True, label="bin5")
    sns.kdeplot(data.loc[data["redshift_group"] == 7]["ABratio"],
                shade=True, label="bin7")
    sns.kdeplot(data.loc[data["redshift_group"] == 9]["ABratio"],
                shade=True, label="bin9")

    plt.xlabel("a/b ratio")
    plt.title("a/b ratio distribution across redshift bins")
    plt.show()
