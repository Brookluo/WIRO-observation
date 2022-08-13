import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval


def show_imstat(line: np.ndarray):
    """Show image statistics. Equivalent to imstat in IRAF.

    Parameters
    ----------
    line : np.ndarray
        the line to show statistics for

    Returns
    -------
    mean, std, median : float
        the mean, standard deviation and median of the image
    """    
    # to prevent overflow for uint16 for adu
    line = line.ravel().astype(np.int64)
    print("Mean:", np.mean(line),
          "Std:", np.std(line),
          "Median:", np.median(line),
          "RMS:", np.sqrt(np.mean(line)),
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
    """A wrapper for plotting an image with zscale.

    Parameters
    ----------
    imdata : np.ndarray
        the image to plot
    ax_handle : matplotlib.axes.Axes
        the axes to plot the image on
    cmap : matplotlib.cmap, optional
        colormap, by default None
    vmin : float, optional
        lower limit for color map, by default None
    vmax : _type_, optional
        upper limit for color map, by default None
    """    
    zscale = ZScaleInterval()
    if not vmin or not vmax:
        vmin, vmax = zscale.get_limits(imdata)
    fig = ax_handle.get_figure()
    im = ax_handle.imshow(imdata, cmap=cmap, vmin=vmin, vmax=vmax, origin='lower')
    if cmap:
        fig.colorbar(im, ax=ax_handle)
    return vmin, vmax
    
    