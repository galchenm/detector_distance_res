import h5py
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
import sys


def get_dataset_keys(f):
    keys = []
    f.visit(lambda key : keys.append(key) if isinstance(f[key], h5py.Dataset) else None)
    return keys

input_path = sys.argv[1] #'/asap3/petra3/gpfs/p11/2021/data/11010507/raw/20210916_6946/'
hdf5path = sys.argv[2] #'entry/data/data'
energy = sys.argv[3] #12 keV
clen = sys.argv[4] #199.6 mm
pixel_size = sys.argv[5] #0.15 mm --> for binned data (binning 2x2) in Eiger case pixel size equals to 150 micron
orgx = sys.argv[6] #1065 - this value you can get from the geometry file
orgy = sys.argv[7] #1092 - this value you can get from the geometry file
extension = sys.argv[8] #'h5'

center = (orgx, orgy)
dpi = 4800

files = glob.glob(os.path.join(input_path, f'*.{extension}'))
radialprofiles = []

for filename in files:
    with h5py.File(filename,'r') as f:
        datasets_names = get_dataset_keys(f)
        print(datasets_names)
        print(hdf5path in datasets_names)
        if hdf5path in datasets_names:
            datasets = f[hdf5path]
            for index in range(0, datasets.shape[0]//100): #data.shape[0]
                data = datasets[index,]
                data = np.where(data>60000,0,data)
                data = np.where(data<0,0,data)          
                y, x = np.indices((data.shape))
                r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
                r = r.astype(int)
                tbin = np.bincount(r.ravel(), data.ravel())
                nr = np.bincount(r.ravel())
                radialprofile = tbin / nr
                radialprofiles.append(radialprofile)


if len(radialprofiles) > 0:
    d = list(map(lambda i: 10 / ((12.4/energy) / (2*np.sin(0.5*np.arctan(i*pixel_size/clen)))), range(1, len(radialprofiles[0])+1)))
    y_est = np.mean(radialprofiles, axis=0)
    y_err = np.std(radialprofiles, axis=0)
    ig, ax = plt.subplots()
    ax.plot(d, y_est, '-', linewidth=1)
    ax.fill_between(d, y_est - y_err, y_est + y_err, linewidth=5, alpha=0.35)
    ax.set_xlabel(r'Q, $nm^{-1}$', labelpad=10)
    ax.set_ylabel(r'I, photons', labelpad=10)
    plt.show()
    #fig.savefig('result.svg', format='svg', dpi=dpi)
