from astropy.io import fits
import os
from Explorer import GalaxyInfo
import bz2
from scipy.interpolate import griddata
import cv2

import numpy as np
# download / discard 'ugriz' file

lib_flag = "/Users/honka/Desktop/galaxy/FITS/"


class getRec:
    """
    Usage: To find each celestial object in the field image
           (256*Ny pixels, 5 bands -ugriz, stored as .bz2 in SDSS),
           it is necessary to get the bound [rmin, rmax, cmin, cmax]
           to conduct cropping from the field image.
           Those bounds are for 2048*1489 pixels image.

           We need to either transform the 256*Ny images into 2048*1489 ones
           by using XINTERP and YINTERP vectors or multiply the bound with
           [256/2048, 256/2048, Ny/1489, Ny/1489].

    Methods:
           1. download the image files by using command line.
              (could be downloaded in any batch size)
           2. unzip.
           3. Get the bound information by querying DBObj: atlasOutline (one by one).
           4. Cropping. Type1: rescale the bound vector;
                        Type2: interpolate the image, back to 2048*1489.

    """
    def __init__(self, galaxy):
        self.galaxy = GalaxyInfo(galaxy=galaxy)
        self.url = self.galaxy.url
        self.file_list = self.getFileName()

    # download by list
    # wget -i /Users/honka/Desktop/galaxy/url.txt
    #  -P /Users/honka/desktop/galaxy

    # file name: /Users/honka/Desktop/galaxy/FITS/
    # frame-r-002662-1-0126.fits.bz2

    def getFITS(self, file):
        decomp = bz2.BZ2File(file)
        hdulist = fits.open(decomp)
        return hdulist

    def getFileName(self):
        list = []
        lib = lib_flag
        for url in self.url:
            list.append(lib + url[-30:])
        return list

    def download(self):
        """
        :return: Download ugriz 5 .bz2 files into lib.
                 No message would be shown, usually would take 30s.
        """
        for url in self.url:
            os.system('wget -q ' + url +
                      ' -P ' + " " + lib_flag)

    def getBound(self):
        """
        :return: Query the database API to get the bound info.
        """
        bound = []
        SQL = "http://skyserver.sdss.org/dr14/en/tools/search/" \
              "x_results.aspx?searchtool=SQL&TaskName=" \
              "Skyserver.Search.SQL&syntax=NoSyntax&ReturnHtml=" \
              "true&cmd=SELECT+objID%2C+rmin%2C+rmax%2C+cmin%2C+cmax%0D%0AFROM" \
              "+atlasOutline%0D%0AWHERE+objID+%3D+"+self.galaxy.galaxy\
              +"&format=html&TableName="

        soup = self.galaxy.query_site(SQL)
        td = soup.body.find_all('td')
        if len(td) != 10:
            raise ValueError("table does not contain required columns")
        else:
            for i in range(6, 10):
                bound.append(float(td[i].text))
        return bound

    def rectangular(self, type='nano'):
        """
        :return: cropped matrix.
                 type1 -- nano
                 type2 -- count
        """
        rectangular = []
        bound = self.getBound()
        cropping = []

        for file in self.file_list:
            text = self.getFITS(file)
            r_bound = [round(bound[0]), round(bound[1])]
            c_bound = [round(bound[2]), round(bound[3])]

            if type == 'nano':
                matrix = text[0].data
                cropping = matrix[r_bound[0]:r_bound[1], c_bound[0]:c_bound[1]]

            if type == 'count':
                data = text[2].data[0]
                xinterp = data[1]
                yinterp = data[2]
                matrix = data[0]
                sky_img = self.back_interpolation(xinterp, yinterp, matrix)

                matrix = text[0].data
                new_matrix = matrix / text[1].data + sky_img
                cropping = new_matrix[r_bound[0]:r_bound[1], c_bound[0]:c_bound[1]]

            rectangular.append(cropping)
        return rectangular

    def back_interpolation(self, xinterp, yinterp, matrix):
        n, m = matrix.shape

        # original grids [n, m]
        x = np.arange(0, n, 1)
        y = np.arange(0, m, 1)
        xx, yy = np.meshgrid(x, y)

        points = np.array((xx.flatten(), yy.flatten())).T  # (49152, 2)
        values = matrix.flatten()  # (49152, )

        xxi, yyi = np.meshgrid(xinterp, yinterp)
        target = np.array((xxi.flatten(), yyi.flatten())).T
        # new = griddata(points, values, (xxi, yyi))
        new = griddata(points, values, target)
        new = new.reshape((2048, 1489))

        return new

    def discard(self):
        os.system("rm -rf /Users/honka/Desktop/galaxy/FITS/")

if __name__ == "__main__":
    rec = getRec(galaxy="1237656906344824834")
    rec.download()
    rect = rec.rectangular()
    print(rect)



