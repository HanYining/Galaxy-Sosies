from astropy.io import fits
import os
from Explorer import GalaxyInfo
import bz2
import scipy.interpolate as interpolate

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

    def rectangular(self, type='type1'):
        """
        :return: cropped matrix.
                 type1 -- rescale bounding vectors
                 type2 -- interpolation image
        """
        rectangular = []
        bound = self.getBound()

        for file in self.file_list:
            text = self.getFITS(file)
            data = text[2].data[0]
            matrix = data[0]
            n, m = matrix.shape

            if type == 'type1':
                # xinterp = data[1]
                # yinterp = data[2]

                r_bound = [round(bound[0]*n/2048), round(bound[1]*n/2048)]
                c_bound = [round(bound[2]*m/1489), round(bound[3]*m/1489)]

                cropping = matrix[r_bound[0]:r_bound[1], c_bound[0]:c_bound[1]]
                rectangular.append(cropping)

            if type == 'type2':
                xinterp = data[1]
                yinterp = data[2]
                new_matrix = self.back_interpolation(xinterp, yinterp, matrix)
        return rectangular

    def back_interpolation(self, xinterp, yinterp, matrix):
        # NOTE: interpolation function needs to be improved
        return interpolate.interp2d(xinterp, yinterp, matrix)

    def discard(self):
        os.system("rm -rf /Users/honka/Desktop/galaxy/FITS/")

if __name__ == "__main__":
    rec = getRec(galaxy="1237656906344824834")
    rec.download()
    rect = rec.rectangular()
    print(rect)



