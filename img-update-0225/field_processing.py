import h5_file
import h5py
import pandas as pd
from getRec import getRec
import numpy as np

# The result comes from SQL Query
galaxy_list = pd.read_csv("SQL_field83.csv")

u = []
g = []
r = []
ii = []
z = []

File_List = []
for i in range(100):
    print("process: " + str(i))

    galaxy = galaxy_list['#Table1'][i+2]
    galaxy_rec = getRec(galaxy=galaxy)

    if galaxy_rec.file_list[0] not in File_List:
        File_List.append(galaxy_rec.file_list[0])
        galaxy_rec.download()

    rect_list = galaxy_rec.rectangular()
    u.append(rect_list[0])
    g.append(rect_list[1])
    r.append(rect_list[2])
    ii.append(rect_list[3])
    z.append(rect_list[4])

# create new h5py file
# don't run this step if there exists one 'galaxy.h5py' file
f = h5py.File("GALAXY.hdf5", 'w')

for i in range(100):
    galaxy = galaxy_list['#Table1'][i + 2]
    ds = f.create_dataset(galaxy, (1,), np.float32)
    ds.attrs['u'] = u[i]
    ds.attrs['g'] = g[i]
    ds.attrs['r'] = r[i]
    ds.attrs['i'] = ii[i]
    ds.attrs['z'] = z[i]

f.close()


