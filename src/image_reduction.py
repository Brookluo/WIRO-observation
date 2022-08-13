import ccdproc as ccdp
from astropy.nddata import CCDData
import astropy.units as u
import numpy as np
import os
from pathlib import Path


def overscan_sub_trim(input_imdata, overscan_fit: callable):
    """Subtract the overscan region and then trim the image to remove
    the overscan region. This is a translation of Chip's WDPzero.cl
    http://physics.uwyo.edu/~WIRO/DoublePrime/WDPzero.cl

    Parameters
    ----------
    input_imdata : numpy.ndarray or CCDData
        the numpy array representing the image
    overscan_fit : callable
        the astropy.modeling.models to fit the overscan region

    Returns
    -------
    output_im : CCDData
        the overscan subtracted and trimmed image
    """
    # Some really important notes about python and IRAF conversions
    # 1. IRAF is 1 based, while python is 0 based
    # 2. IRAF [x:y] is inclusive on both ends, ie, it will have (y-x+1) items
    #    while python [x:y) is inclusive on left, not on tight, thus have (y-x) items
    #    in IRAF img[1:200] == in python img[0:200] ([left -1, right -1 +1], shift left
    #    and then +1 to include right element
    # 3. IRAF image is the transpose (flip) of the python image
    # Dimensions in fits header = (column, row)
    # shape of numpy array = (row, column)
    if not isinstance(input_imdata, CCDData):
        if isinstance(input_imdata, np.ndarray):
            ccd_im = CCDData(input_imdata, unit=u.adu)
        else:
            ccd_im = CCDData.read(input_imdata, unit=u.adu)
    else:
        # to avoid any potential inplace operation contamination
        ccd_im = input_imdata.copy()
    # ! needs to pay extra attention to the indexing method (xy or ij), column-row or row-column
    # Now it is row-column, transpose the image to get column-row
    ccd_im.data = ccd_im.data.T
    # Following WDPzero.cl for four amplifiers
    amp1 = ccd_im[8 - 1 : 2100, 0 : 2048]
    amp2 = ccd_im[2101 - 1 : 4193, 0 : 2048]
    amp3 = ccd_im[8 - 1 : 2100, 2049 - 1 : 4096]
    amp4 = ccd_im[2101 - 1 : 4193, 2049 - 1 : 4096]
    amp1.data = amp1.data.T
    amp2.data = amp2.data.T
    amp3.data = amp3.data.T
    amp4.data = amp4.data.T
    
    # overscan can only subtract columns, so need to use all rows
    amp1_z = ccdp.subtract_overscan(
        amp1, overscan=amp1[:, 2054 - 1 : 2089], model=overscan_fit
    )
    amp2_z = ccdp.subtract_overscan(
        amp2, overscan=amp2[:, 4 - 1 : 40], model=overscan_fit
    )
    amp3_z = ccdp.subtract_overscan(
        amp3, overscan=amp3[:, 2054 - 1 : 2088], model=overscan_fit
    )
    amp4_z = ccdp.subtract_overscan(
        amp4, overscan=amp4[:, 4 - 1 : 40], model=overscan_fit
    )
    amp1_z.data = amp1_z.data.T
    amp2_z.data = amp2_z.data.T
    amp3_z.data = amp3_z.data.T
    amp4_z.data = amp4_z.data.T

    amp1_zt = ccdp.trim_image(amp1_z[0 : 2048, :])
    amp2_zt = ccdp.trim_image(amp2_z[46 - 1 : 2093, :])
    amp3_zt = ccdp.trim_image(amp3_z[0 : 2048, :])
    amp4_zt = ccdp.trim_image(amp4_z[46 - 1 : 2093, :])

    # to preserve the output shape
    output_im = CCDData(np.ones((4096, 4096)), unit=u.adu)
    output_im.data[0 : 2048, 0 : 2048] = amp1_zt.data
    output_im.data[2049 - 1 : 4096, 0 : 2048] = amp2_zt.data
    output_im.data[0 : 2048, 2049 - 1 : 4096] = amp3_zt.data
    output_im.data[2049 - 1 : 4096, 2049 - 1 : 4096] = amp4_zt.data
    
    output_im.data = output_im.data.T
    
    return output_im


def inv_median(a):
    return 1 / np.median(a)


def make_masterflat(filelist, output_dir, band, overwrite=False):
    """Make a master flat for a given band using median combine with 3-sigma clipping.
    Two master flats are created. One with normalization by mean of the master flat, one
    is without the normalization.

    Parameters
    ----------
    filelist : list
        a list of filenames of the given band flat images. All filenames must be
        full path to that file.
    output_dir : str
        the output directory to save the master flat files.
    band : str
        a string representing the band. (UBVRI or ugriz)
    overwrite : bool, optional
        if a masterflat already exists, whether overwrite that file, by default False
        
    Returns
    -------
    CCDdata
        a CCDData object of the normalized, count weighted master flat
    """
    from utils import plot_zscale_image, show_imstat
    import matplotlib.pyplot as plt

    # use inverse median to scale all the image to unity first
    #         flat_d.data /= np.median(flat_v_d.data)
    #         flat_d.write(f"{root_dir}/a{i:0>3}_d_medscaled.fits", overwrite=True)
    # use original image brightness
    mean_count = np.array([np.mean(CCDData.read(file).data) for file in filelist])
    mean_count /= np.sum(mean_count)
    combined_flat_clip_med_weighted_avg = ccdp.combine(
        filelist,
        method="median",
        weights=mean_count,
        scale=inv_median,
        sigma_clip=True,
        sigma_clip_low_thresh=3,
        sigma_clip_high_thresh=3,
        sigma_clip_func=np.ma.mean,
        sigma_clip_dev_func=np.ma.std,
    )
    combined_flat_clip_med_weighted_avg.write(
        os.path.join(output_dir, f"masterflat{band}_clip_med_weighted_count.fits"),
        overwrite=overwrite,
    )
    combined_flat_clip_med_weighted_avg.data /= np.mean(
        combined_flat_clip_med_weighted_avg.data
    )
    combined_flat_clip_med_weighted_avg.write(
        os.path.join(output_dir, f"masterflat{band}_norm.fits"), overwrite=overwrite
    )
    fig, ax = plt.subplots(1, 2, figsize=(30, 10))
    plot_zscale_image(combined_flat_clip_med_weighted_avg.data, ax[0], "gray")
    ax[0].set_aspect("equal")
    #     ax.set_title(f"{masterVs_name[i]}")
    ax[1].plot(combined_flat_clip_med_weighted_avg.data[:, 200].ravel())
    ax[1].set_ylabel("counts")
    ax[1].set_xlabel("pixel")
    show_imstat(combined_flat_clip_med_weighted_avg.data)
    return combined_flat_clip_med_weighted_avg 


def reduce_images(
    sci_img: list,
    bias_img: list,
    dark_img: list,
    flat_img: dict[str, list],
    output_dir: str,
    overwrite=False,
):
    """Reduce a list of images using the master flat and master bias.
    
    """
    from astropy.modeling import models, fitting, polynomial

    # first do overscan subtraction on all images
    all_img = sci_img
    if bias_img:
        all_img.extend(bias_img)
    if dark_img:
        all_img.extend(dark_img)
    if flat_img:
        for v in flat_img.values():
            all_img.extend(v)
    # fit_poly = fitting.LevMarLSQFitter()
    poly = polynomial.Polynomial1D(degree=3)
    cheb = polynomial.Chebyshev1D(degree=3)
    leg = polynomial.Legendre1D(degree=3)
    herm = polynomial.Hermite1D(degree=3)
    
    all_img_zt = []
    for img in all_img:
        ccd_im = CCDData.read(img, unit=u.adu)
        # Default to use a chebyshev polynomial
        ccd_zt = overscan_sub_trim(ccd_im, overscan_fit=cheb)
        file_loc = Path(output_dir, img.stem + "_zt.fits")
        ccd_zt.write(file_loc, overwrite=overwrite)
        all_img_zt.append(file_loc)
    if bias_img:
        # make master bias and bias subtraction
        bias_zt_files = [os.path.join(output_dir, img.stem + "_zt.fits") for img in bias_img]
        combined_bias = ccdp.combine(
            bias_zt_files,
            method="average",
            sigma_clip=True,
            sigma_clip_low_thresh=5,
            sigma_clip_high_thresh=5,
            sigma_clip_func=np.ma.median,
            sigma_clip_dev_func=np.ma.std,
            mem_limit=350e6,
        )
        combined_bias.write(os.path.join(output_dir, "masterbias.fits"), overwrite=overwrite)
    # TODO need a criterion to decide whether to subtract masterbias or not
    for img in all_img_zt:
        ccd_im = CCDData.read(img, unit=u.adu)
        ccd_b = ccdp.subtract_bias(ccd_im, combined_bias)
        # from _zt to _b suffix
        ccd_b.write(os.path.join(output_dir, img.stem.rsplit("_")[0] + "_b.fits"), overwrite=overwrite)
    
    if dark_img:
        # make master dark
        dark_b_files = [img.stem + "_zt.fits" for img in dark_img]
        # TODO check the exposure time, use the longest exposure ones
        combined_dark_clip_med = ccdp.combine(
            dark_b_files,
            method='median',
            sigma_clip=True, sigma_clip_low_thresh=3, sigma_clip_high_thresh=3,
            sigma_clip_func=np.ma.mean, sigma_clip_dev_func=np.ma.std,
            mem_limit=350e6
        )
        combined_dark_clip_med.write(os.path.join(output_dir, "masterdark.fits"), overwrite=overwrite)
        # For WIRO, dark current is not a significant issue, so we need to think
        # whether remove the dark current or not.
    if flat_img:
        for band, files in flat_img.items():
            make_masterflat(files, output_dir, band, overwrite=overwrite)
    # TODO Distortion correction
    # linearity check