import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
import ccdproc as ccdp
import astropy.units as u


def show_imstat(line: np.ndarray):
    line = line.ravel()
    print("Mean:", np.mean(line),
          "Std:", np.std(line),
          "Median:", np.median(line),
          "RMS:", np.sqrt(np.mean(line ** 2)),
          "Min:", np.min(line),
          "Max:", np.max(line)
          )
    return np.mean(line), np.std(line), np.median(line)


def plot_line(line: np.ndarray):
    plt.figure(figsize=(8, 6))
    plt.plot(line.ravel())
    plt.ylabel("counts")
    plt.xlabel("pixel")
    plt.show()
    
    
def plot_zscale_image(imdata: np.ndarray, ax_handle, cmap=None, vmin=None, vmax=None):
    zscale = ZScaleInterval()
    if not vmin or not vmax:
        vmin, vmax = zscale.get_limits(imdata)
    fig = ax_handle.get_figure()
    im = ax_handle.imshow(imdata, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')
    if cmap:
        fig.colorbar(im, ax=ax_handle)
    
    
def overscan_subs_trim(imdata, overscan_fit):
    if type(imdata) is not ccdp.CCDData:
        ccd_im = ccdp.CCDData(imdata, unit=u.adu)
    # these 50 columns should be good representatives for the overscan region
    ccd_z = ccdp.subtract_overscan(ccd_im, overscan=ccd_im[:, 2120:2170], model=overscan_fit)
    # science image should be good in columns 55:2100
    ccd_zt = ccdp.trim_image(ccd_z[:, 55:2100])
    return ccd_zt