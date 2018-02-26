import bs4
import requests

# explore a galaxy by u
# sing its ID


class GalaxyInfo(object):
    """
    Usage: Get u, g, r, i, z URLs for one celestial object.
    Method:
            SDSS provides a simple method called 'dbo.fGetUrlSpecImg()'.
            Some objects can not be found by using this method
            (but more often we can make a query like this).

            This function is written to query the navigator and get the URLs for
            the fits files of the five bands.

            Searching by using navigator:
            1. go to http://skyserver.sdss.org/dr14/en/tools/chart/navi.aspx.
            2. click 'explore'.
            3. click 'FITS'.
            4. click 'Corrected Frames' to see the description of the FITS file.
            As described, field images are saved as ALLSKY	float32[256,ny].
            Interpolation vectors (bilinear interp) are stored as XINTERP and YINTERP.

    """
    def __init__(self, galaxy):
        """
        :param galaxy: galaxy id - string
        :attrs: site - explore site in step 2.
                profile - site's content.
                fits_site - FITS website in step 3.
                fits - FITS site's content
                ugriz - we want files about all the five ugriz bands.
                url - URL links to the database.
        """
        self.galaxy = galaxy
        self.site = "http://skyserver.sdss.org/dr14/en/tools/" \
                    "explore/Summary.aspx?id=" + galaxy
        self.profile = self.query_site(self.site)
        self.content = self.profile.body.find(self._get_FitsWeb)

        self.fits_site = "http://skyserver.sdss.org/dr14/en" \
                         "/tools/explore/" + self.content['href']

        self.fits = self.query_site(self.fits_site)
        self.ugriz = ['u', 'g', 'r', 'i', 'z']
        self.url = self.get_ugriz()

    def query_site(self, site):
        """
        :param site: website url
        :return: soup object for site
        """
        try:
            web = requests.get(site)
        except requests.exceptions.RequestException as e:
            raise ValueError('page not found')

        soup = bs4.BeautifulSoup(web.text, 'lxml')
        return soup

    def get_FitsWeb(self):
        """
        :return: get url to the FITS summary website
        """
        return self.profile.body.find(self._get_FitsWeb)

    def _get_FitsWeb(self, tag):
        return tag.has_attr('title') and "Get FITS images " \
                                         "of the SDSS fields containing this object." \
                                         in tag['title']

    def _get_ugriz(self, key):
        """
        Help to return all the links for corrected frames,
        binned frames, and masked frames, in case they are needed
        in the future.
        """
        table_content = self.fits.find_all('td')
        index_list = []

        for i, element in enumerate(table_content):
            content = element.get_text(strip=True)
            if key == content:
                index_list.append(i)
        if len(index_list) > 0:
            return index_list  # [corrected, binned, masked]
        else:
            return -2

    def get_ugriz(self):
        """
        :return: 5 URLs for ugriz file in the database.
                 Files are stored as .bz2 (zipped).
        """
        fits_url = []
        for key in self.ugriz:
            index = self._get_ugriz(key)
            if index == -2:
                raise ValueError("no fits index_list")
            else:
                # only return the corrected frames
                fits_url.append(self.fits.find_all('td')
                                [index[0]].find('a').attrs['href'])
        return fits_url


if __name__ == "__main__":
    target = GalaxyInfo(galaxy="1237656906344824833")
    url = target.url
    print(url)




