from astropy.io import fits
import numpy as np
import glob
import random
import gc
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pyplot as plt
from scipy import stats


class galaxy_processor():
    '''
    A galaxy processor to filter through the fits file and get galaxies.
    '''

    def __init__(self, loc="data/*.fits"):
        self.file_names = glob.glob(loc)
        self.galaxy_cnt = 0
        self.cnt = 0
        self.flux, self.wav = self._filter_galaxy()
        self.flux, threshold = self._align_spec(self.flux, self.wav)
        self.loglam_min, self.loglam_max = threshold
        self.wavelen = self._generate_wavelength()

    def _filter_galaxy(self):
        res = []
        wav = []
        for file in self.file_names:
            if self.cnt > 0 and self.cnt % 50 == 0:
                print("Successfully Processed: " + str(self.cnt))
                print("Total galaxy: " + str(self.galaxy_cnt))

            with fits.open(file, memmap=False) as fits_file:
                if fits_file[2].data["class"][0] == "GALAXY":
                    res.append(fits_file[1].data["model"])
                    wav.append(fits_file[1].data["loglam"])
                    self.galaxy_cnt += 1
                self.cnt += 1
                gc.collect()
        return res, wav

    def _align_spec(self, res, wav):
        maxx = min([wavk[-1] for wavk in wav])
        minn = max([wavk[0] for wavk in wav])
        res_new = []
        for i, resk in enumerate(res):
            res_new.append(res[i][(wav[i] > minn) & (wav[i] < maxx)])
        return np.vstack(res_new), [minn, maxx]

    def _generate_wavelength(self):
        return 10 ** np.linspace(self.loglam_min,
                                 self.loglam_max,
                                 self.flux.shape[1])

    def extract_visible(self, start, end):
        ind = (self.wavelen >= start) & (self.wavelen <= end)
        return self.flux[:, ind], self.wavelen[ind]

    def detect_balmer(self, local_region=20):
        balmer_wavelen = 3650
        left = balmer_wavelen - local_region
        right = balmer_wavelen + local_region
        left_side, waves = self.extract_visible(left, balmer_wavelen)
        right_side, waves = self.extract_visible(balmer_wavelen, right)
        self.balmer_diff = np.mean(left_side, axis=1) - \
            np.mean(right_side, axis=1)

        return self.balmer_diff

    def detect_lentune(self):

        fig = plt.figure()
        density = stats.kde.gaussian_kde(self.balmer_diff)
        x = np.linspace(-0.05, 0.05, 2000)
        density.covariance_factor = lambda: .01
        density._compute_covariance()
        plt.plot(x, density(x))
        fig.show()

        pass


def _plot(wave, x, mph, mpd, threshold, edge, valley, ax, ind):
    """Plot results of the detect_peaks function, see its help."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib is not available.')
    else:
        if ax is None:
            _, ax = plt.subplots(1, 1, figsize=(8, 4))

        ax.plot(wave, x, 'b', lw=1)
        if ind.size:
            label = 'valley' if valley else 'peak'
            label = label + 's' if ind.size > 1 else label
            ax.plot(wave[ind], x[ind], '+', mfc=None, mec='r', mew=2, ms=8,
                    label='%d %s' % (ind.size, label))
            ax.legend(loc='best', framealpha=.5, numpoints=1)
        ymin, ymax = x[np.isfinite(x)].min(), x[np.isfinite(x)].max()
        yrange = ymax - ymin if ymax > ymin else 1
        ax.set_ylim(ymin - 0.1*yrange, ymax + 0.1*yrange)
        ax.set_xlabel('wavelength', fontsize=14)
        ax.set_ylabel('flux', fontsize=14)
        mode = 'Valley detection' if valley else 'Peak detection'
        ax.set_title("%s (mph=%s, mpd=%d, threshold=%s, edge='%s')"
                     % (mode, str(mph), mpd, str(threshold), edge))
        # plt.grid()
        plt.show()


def detect_peaks(wave, x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    """
    Detect peaks in data based on their amplitude and other features.
    Parameters
    x : 1D array
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns:
    ind : 1D array_like indeces of the peaks in `x`.
    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    if show:
        if indnan.size:
            x[indnan] = np.nan
        if valley:
            x = -x
        _plot(wave, x, mph, mpd, threshold, edge, valley, ax, ind)

    return ind


def randomplot(flux, wave, num, lower, upper):

    fig, axes = plt.subplots(num, 1)
    random_indx = [random.randint(0, flux.shape[0])
                   for _ in range(num)]

    for i, ind in enumerate(random_indx):
        axes[i].plot(wave, flux[ind, ])
        axes[i].set_ylim([lower, upper])

    fig.show()
    input("<Hit Enter To Close>")
    plt.close(fig)
    pass


test = galaxy_processor()
test1 = test.detect_balmer(10)

plt.plot(test.wavelen, test.flux[15, ])
plt.title("a clear Plateau after 600 nm")
plt.xlabel("wave length")
plt.ylabel("flux")
plt.show()


plt.plot(test.wavelen, test.flux[35, ])
plt.title("Spectrum without a Plateau after 600 nm")
plt.xlabel("wave length")
plt.ylabel("flux")
plt.show()

fig = plt.figure()
density = stats.kde.gaussian_kde(test.balmer_diff)
x = np.linspace(-0.05, 0.05, 2000)
density.covariance_factor = lambda: .01
density._compute_covariance()
plt.plot(x, 100*density(x)/sum(density(x)))
plt.xlabel("difference")
plt.ylabel("galaxy %")
plt.title("Kernel smoothed distribution plots")
fig.show()


peaks = detect_peaks(test.flux[1, ], np.mean(test.flux[1, ]))


sub_flux, sub_wave = test.extract_visible(3600, 3800)
randomplot(test.flux, test.wavelen, 2, 0, 5)

plt.plot(sub_wave, sub_flux[random.randint(0, sub_flux.shape[0]-1), ])
plt.show()

plt.plot(test.wavelen, test.flux[2, ])
plt.show()

# here we do a PCA decomposition to see how the top trend of the
# data looks like.
spec_norm = preprocessing.scale(spec)
plt.plot(spec_norm[0,])
plt.show()
pca = PCA()
pca.fit(spec_norm)
first_component = pca.components_[0, ]
second_component = pca.components_[1, ]
third_component = pca.components_[2, ]


plt.plot(first_component, "r")
plt.show()


test_file = fits.open("/Users/yinhan/stats_prac/data/spec-3586-55181-0104.fits")
hdu1, hdu2, hdu3, hdu4 = test_file
test_file2 = fits.open("/Users/yinhan/stats_prac/data/spec-3586-55181-0302.fits")
hdu21, hdu22, hdu23, hdu24 = test_file2
