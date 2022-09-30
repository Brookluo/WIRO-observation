import ccdproc as ccdp
from astropy.nddata import CCDData
import astropy.units as u
import numpy as np
from pathlib import Path
import argparse
import sys

        
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
    # type(adu) == uint16, what if we have minus counts causing overflow? (i.e.)
    # when the image count is actually all noise
    # change the datatype to int 64 to prevent that
    ccd_im.data = ccd_im.data.astype(np.int64)
    # Now it is row-column, transpose the image to get column-row
    ccd_im.data = ccd_im.data.T
    # Following WDPzero.cl for four amplifiers
    amp1 = ccd_im[8 - 1 : 2100, 0:2048]
    amp2 = ccd_im[2101 - 1 : 4193, 0:2048]
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

    amp1_zt = ccdp.trim_image(amp1_z[0:2048, :])
    amp2_zt = ccdp.trim_image(amp2_z[46 - 1 : 2093, :])
    amp3_zt = ccdp.trim_image(amp3_z[0:2048, :])
    amp4_zt = ccdp.trim_image(amp4_z[46 - 1 : 2093, :])

    # to preserve the output shape
    output_im = CCDData(np.ones((4096, 4096)), unit=u.adu)
    output_im.data[0:2048, 0:2048] = amp1_zt.data
    output_im.data[2049 - 1 : 4096, 0:2048] = amp2_zt.data
    output_im.data[0:2048, 2049 - 1 : 4096] = amp3_zt.data
    output_im.data[2049 - 1 : 4096, 2049 - 1 : 4096] = amp4_zt.data

    output_im.data = output_im.data.T
    output_im.data = output_im.data.astype(np.int16)
    output_im.header = ccd_im.header

    return output_im


def all_overscan_sub_trim(filelist: list, output_dir: str, polyfit='cheb', 
                          overwrite=False):
    """_summary_

    Parameters
    ----------
    filelist : list
        _description_
    output_dir : str
        _description_
    polyfit : str, optional
        _description_, by default 'cheb'
    overwrite : bool, optional
        _description_, by default False

    Raises
    ------
    ValueError
        _description_
    """    
    from astropy.modeling import polynomial
    
    # list a bunch of polynomial models, usually 3rd order is enough
    poly = polynomial.Polynomial1D(degree=3)
    cheb = polynomial.Chebyshev1D(degree=3)
    leg = polynomial.Legendre1D(degree=3)
    herm = polynomial.Hermite1D(degree=3)
    if polyfit == 'cheb':
        fit_func = cheb
    elif polyfit == 'poly':
        fit_func = poly
    elif polyfit == 'leg':
        fit_func = leg
    elif polyfit == 'herm':
        fit_func = herm
    else:
        raise ValueError(f"Unknown polynomial overscan fitting model {polyfit}")
        
    # assume all images have the same suffix
    file_suffix = filelist[0].suffix
    cur_stage = "z"
    for img in filelist:
        ccd_im = CCDData.read(img, unit=u.adu)
        # Default to use a chebyshev polynomial
        ccd_zt = overscan_sub_trim(ccd_im, overscan_fit=fit_func)
        post_zt_file_loc = Path(output_dir, f"{img.stem}_{cur_stage}{file_suffix}")
        ccd_zt.write(post_zt_file_loc, overwrite=overwrite)


def bias_subtract(filelist: list, bias_filelist:list, output_dir: str, 
                  overwrite=False):
    """_summary_

    Parameters
    ----------
    filelist : list
        _description_
    bias_filelist : list
        _description_
    output_dir : str
        _description_
    overwrite : bool, optional
        _description_, by default False
    """
    # assume all images have the same suffix
    file_suffix = filelist[0].suffix
    # make master bias and bias subtraction
    # bias_zt_files = [Path(output_dir, f"{img.stem}_{cur_stage}{file_suffix}") for img in bias_filelist]
    combined_bias = ccdp.combine(
        bias_filelist,
        method="average",
        sigma_clip=True,
        sigma_clip_low_thresh=5,
        sigma_clip_high_thresh=5,
        sigma_clip_func=np.ma.median,
        sigma_clip_dev_func=np.ma.std,
        mem_limit=350e6,
        unit=u.adu
    )
    combined_bias.write(Path(output_dir, "masterbias.fits"), overwrite=overwrite)
    # TODO need a criterion to decide whether to subtract masterbias or not
    
    cur_stage = "b"
    for img in filelist:
        if img in bias_filelist:
            continue
        ccd_im = CCDData.read(img, unit=u.adu)
        ccd_b = ccdp.subtract_bias(ccd_im, combined_bias)
        # from _z to _zb suffix
        if "_" in img.stem:
            new_fname = f"{img.stem+cur_stage}{file_suffix}"
        else:
            new_fname = f"{img.stem}_{cur_stage}{file_suffix}"
        post_bias_file_loc = Path(output_dir, new_fname )
        ccd_b.write(post_bias_file_loc, overwrite=overwrite)


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
    from plot_utils import plot_zscale_image, show_imstat
    import matplotlib.pyplot as plt

    # use inverse median to scale all the image to unity first
    #         flat_d.data /= np.median(flat_v_d.data)
    #         flat_d.write(f"{root_dir}/a{i:0>3}_d_medscaled.fits", overwrite=True)
    # use original image brightness
    mean_count = np.array([np.mean(CCDData.read(file, unit=u.adu).data) for file in filelist])
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
        unit=u.adu
    )
    combined_flat_clip_med_weighted_avg.write(
        Path(output_dir, f"masterflat{band}_clip_med_weighted_count.fits"),
        overwrite=overwrite,
    )
    combined_flat_clip_med_weighted_avg.data /= np.mean(
        combined_flat_clip_med_weighted_avg.data
    )
    combined_flat_clip_med_weighted_avg.write(
        Path(output_dir, f"masterflat{band}_norm.fits"), overwrite=overwrite
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


def flat_correct(filelist: list, masterflat_filelist: dict[str, list], 
                output_dir: str, overwrite=False):
    """_summary_

    Parameters
    ----------
    filelist : list
        _description_
    masterflat_filelist : dict[str, list]
        _description_
    output_dir : str
        _description_
    overwrite : bool, optional
        _description_, by default False
    """    
    import fitsio
    
    file_suffix = filelist[0].suffix
    # separate files into different filters
    filter_image = {band: [] for band in masterflat_filelist.keys()}
    for img in filelist:
        header = fitsio.read_header(img)
        # filter name is the string after column
        # Filter 5: i' 54605
        fitlername = header["FILTER"].rsplit(":")[-1].strip()
        filter_image[fitlername].append(img)
    cur_stage = "f"
    for band, masterflat in masterflat_filelist.items():
        for img in filter_image[band]:
            ccd_im = CCDData.read(img, unit=u.adu)
            ccd_masterflat = CCDData.read(masterflat, unit=u.adu)
            ccd_f = ccdp.flat_correct(ccd_im, ccd_masterflat)
            if "_" in img.stem:
                new_fname = f"{img.stem+cur_stage}{file_suffix}"
            else:
                new_fname = f"{img.stem}_{cur_stage}{file_suffix}"
            post_flat_file_loc = Path(output_dir, new_fname)
            ccd_f.write(post_flat_file_loc, overwrite=overwrite)

# really should decompose this function into smaller functions
# so each function just process one stage
def reduce_images(
    sci_img: dict[str, list],
    bias_img: list,
    dark_img: list,
    flat_img: dict[str, list],
    output_dir: str,
    overwrite=False,
):
    """Reduce a list of images using the master flat and master bias.
    All paths in the arguments must be full path to files. The path 
    should be python Path object.

    Parameters
    ----------
    sci_img : dict[str, list[Path]]
        A dictionary of path to science images. The key is the filter name.
        the value is the list of path to science images for that band.
    bias_img : list[Path]
        A list of path to bias images.
    dark_img : list[Path]
        A list of path to dark images.
    flat_img : dict[str, list[Path]]
        A dictionary of path to flat images. The key is the filter name.
        the value is the list of path to flat images for that filter.
    output_dir : str
        The output directory to save the reduced images.
    overwrite : bool, optional
        whether to overwrite the existing files if any, by default False
    """
    from astropy.modeling import polynomial

    # first do overscan subtraction on all images
    all_img = list(sci_img.values())
    if bias_img:
        all_img.extend(bias_img)
    if dark_img:
        all_img.extend(dark_img)
    if flat_img:
        all_img.extend(list(flat_img.values()))
    poly = polynomial.Polynomial1D(degree=3)
    cheb = polynomial.Chebyshev1D(degree=3)
    leg = polynomial.Legendre1D(degree=3)
    herm = polynomial.Hermite1D(degree=3)
    # assume all images have the same suffix
    file_suffix = all_img[0].suffix
    
    cur_stage = "z"
    for img in all_img:
        ccd_im = CCDData.read(img, unit=u.adu)
        # Default to use a chebyshev polynomial
        ccd_zt = overscan_sub_trim(ccd_im, overscan_fit=cheb)
        post_zt_file_loc = Path(output_dir, f"{img.stem}_{cur_stage}{file_suffix}")
        ccd_zt.write(post_zt_file_loc, overwrite=overwrite)
        
    if bias_img:
        # make master bias and bias subtraction
        bias_zt_files = [Path(output_dir, f"{img.stem}_{cur_stage}{file_suffix}") for img in bias_img]
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
        combined_bias.write(Path(output_dir, "masterbias.fits"), overwrite=overwrite)
        # TODO need a criterion to decide whether to subtract masterbias or not
        last_stage = cur_stage
        cur_stage += "b"
        for img in all_img:
            if img in bias_img:
                continue
            img_last_stage = Path(output_dir, f"{img.stem}_{last_stage}{file_suffix}")
            ccd_im = CCDData.read(img_last_stage, unit=u.adu)
            ccd_b = ccdp.subtract_bias(ccd_im, combined_bias)
            # from _zt to _b suffix
            post_bias_file_loc = Path(output_dir, f"{img.stem}_{cur_stage}{file_suffix}")
            ccd_b.write(post_bias_file_loc, overwrite=overwrite)
            all_img.append(post_bias_file_loc)
            
    if dark_img:
        # make master dark
        dark_b_files = [Path(output_dir, f"{img.stem}_{cur_stage}{file_suffix}") for img in dark_img]
        # TODO check the exposure time, use the longest exposure ones
        combined_dark_clip_med = ccdp.combine(
            dark_b_files,
            method="median",
            sigma_clip=True,
            sigma_clip_low_thresh=3,
            sigma_clip_high_thresh=3,
            sigma_clip_func=np.ma.mean,
            sigma_clip_dev_func=np.ma.std,
            mem_limit=350e6,
        )
        combined_dark_clip_med.write(Path(output_dir, "masterdark.fits"), overwrite=overwrite)
        # For WIRO, dark current is not a significant issue, so we need to think
        # whether remove the dark current or not.
        last_stage = cur_stage
        cur_stage += "d"
        for img in all_img:
            if img in dark_img:
                continue
            img_last_stage = Path(output_dir, f"{img.stem}_{last_stage}{file_suffix}")
            ccd_im = CCDData.read(img_last_stage, unit=u.adu)
            ccd_d = ccdp.subtract_dark(ccd_im, combined_dark_clip_med)
            # from _zt to _b suffix
            dark_file_loc = Path(output_dir, f"{img.stem}_{cur_stage}{file_suffix}")
            ccd_d.write(dark_file_loc, overwrite=overwrite)

    if flat_img:
        for band, files in flat_img.items():
            make_masterflat(files, output_dir, band, overwrite=overwrite)
        # TODO flat field correction on all images
    # TODO Distortion correction
    # linearity check

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Give a observation log file to process the images.")
    parser.add_argument("obs_log_path", help="The observation log file.")
    args = parser.parse_args()
    # obs_log.read(args.obs_log_path)
    
    
