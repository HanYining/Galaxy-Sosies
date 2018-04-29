import h5py
import numpy as np
from getRec import getRec


def h5_galaxy_in(rect_list, galaxy):
    """
    Usage: store the cropping image into h5py file.
    Format: h5py file name 'GALAXY'.
            'GALAXY' contains sub datasets,
            which are named by the galaxy id.
            Each dataset contains 5 columns,
            cropping images in 5 bands are stored inside.
    """
    h5_file = h5py.File('GALAXY', 'r+')
    ds = h5_file.create_dataset(galaxy, (1,), np.float32)
    try:
        ds.attrs['u'] = rect_list[0]
        ds.attrs['g'] = rect_list[1]
        ds.attrs['r'] = rect_list[2]
        ds.attrs['i'] = rect_list[3]
        ds.attrs['z'] = rect_list[4]

        h5_file.close()

    except Exception:
        print("error at " + galaxy)
        h5_file.close()


if __name__ == "__main__":
    galaxy = "1237659930535396672"

    galaxy_rec = getRec(galaxy=galaxy)
    galaxy_rec.download()
    rect_list = galaxy_rec.rectangular()
    h5_galaxy_in(rect_list, galaxy)

    # read the stored image
    # h5_file = h5py.File("GALAXY", 'r')
    # h5_file["/1237659930535396672"].attrs['u']
    # h5_file.close()

    galaxy_rec.discard()



