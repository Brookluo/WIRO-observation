import numpy as np


def convert_flux_to_mag(flux):
    """Convert flux in nanomaggies to magnitude
    follow the transformation in the following website
    https://www.legacysurvey.org/dr9/description/#photometry
    
    Args:
        flux (float): flux of the target in nanomaggies
    Returns:
        float: magnitude of the target
    """    
    m = 22.5 - 2.5 * np.log10(flux)
    return m


def find_sweep_with_ra_dec(ra, dec):
    '''find the sweep file given an ra and dec
    
    Args:
        ra (float): right acesion in degrees
        dec (float): declination in degrees
        
    Return:
        (string): a file name contains this ra and dec 
    '''
    ra_floor = int(ra // 10 * 10)
    dec_floor = int(dec // 5 * 5)
    ra_roof = int(ra_floor + 10)
    dec_roof = int(dec_floor + 5)
    assign_pm = lambda num: 'm' if num < 0 else 'p'
#     floor_prefix = 'p'
#     roof_prefix = 'p'
#     if dec_floor < 0:
#         floor_prefix = 'm'
#     if dec_roof < 0:
#         roof_prefix = 'm'
    floor_prefix = assign_pm(dec_floor)
    roof_prefix = assign_pm(dec_roof)
    filename = f"sweep-{ra_floor:03}{floor_prefix}{abs(dec_floor):03}-{ra_roof:03}{roof_prefix}{abs(dec_roof):03}.fits"
    return filename
