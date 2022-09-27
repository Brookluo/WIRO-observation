import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from astropy.units.quantity import Quantity
from astroplan.scheduling import Scheduler
import numpy as np


def output_minsecz_in_month(month: int, year: int, obs_loc: EarthLocation, obj_coords: SkyCoord):
    """A function used to find objects with minimum airmass (secz)
    of each night in a given month.

    Parameters
    ----------
    month : int
        Month of the year. (1-12)
    year : int
        Year
    obs_loc : EarthLocation
        Location of the observatory.
    obj_coords : SkyCoord
        Single or a list of SkyCoords of the objects to be observed.

    Returns
    -------
    date_tab : Table
        Table with minimum airmass (secz) of each night in a given month.
    """
    import pandas as pd
    import numpy as np
    
    
    # YL setting up time and offset to UTC
    date_dict = []
    # set meridian
    t11 = Time(f"{year}-{month:>02}-01T00:00:00", format="isot", scale="utc")
    utcoffset = -7 * u.hour  # MST
    while t11.ymdhms["month"] == month:
        utc11 = t11 - utcoffset
        field_altaz = obj_coords.transform_to(AltAz(obstime=utc11, location=obs_loc))
        if not obj_coords.shape:
            date_dict.append(
                (
                    t11.value,
                    obj_coords,
                    round(obj_coords.ra.degree, 6),
                    round(obj_coords.dec.degree, 6),
                    round(field_altaz.secz.value, 6),
                )
            )
        else:
            pos_airmass = field_altaz.secz > 0
            idx = np.argmin(field_altaz.secz[pos_airmass])
            date_dict.append(
                (
                    t11.value,
                    obj_coords[pos_airmass][idx][0],
                    round(obj_coords.ra.degree[pos_airmass][idx], 6),
                    round(obj_coords.dec.degree[pos_airmass][idx], 6),
                    round(field_altaz.secz[pos_airmass][idx].value, 6),
                )
            )
        t11 += 1 * u.day
    # YL IO part, write to csv files
    date_tab = pd.DataFrame(
        date_dict,
        columns=[
            "Date",
            "Object Coordinates (hms.ss degree)",
            "RA(o)",
            "Dec(o)",
            "Airmass",
        ],
    )
    # date_tab.to_csv(f"mon_{month}_qso_min_airmass_month.csv", index=False)
    return date_tab


def schedule_obs_for_one_day(month: int, day: int, targets: SkyCoord, scheduler: Scheduler, utcoffset: Quantity, exp: Quantity,
                             n_exposures: int, read_out: Quantity, priority_filters: list[tuple[int, str]]):
    """Schedule observations for one day.

    Parameters
    ----------
    month : int
        Month to schedule observations for. (1-12)
    day : int
        Day to schedule observations for.
    targets : SkyCoord
        SkyCoords of targets to observe.
    scheduler : Scheduler
        Scheduler to use. (PriorityScheduler or SequentialScheduler)
    utcoffset : Quantity
        UTC offset of the observatory.
    exp : Quantity
        Exposure time.
    n_exposures : int
        Number of exposures to take.
    read_out : Quantity
        Readout time.
    priority_filters : list[tuple[int, str]]
        List of Tuples (priority, filters) to use. Recommend using
        list(enumerate(filters)) to generate this.

    Returns
    -------
    day_schedule : Table
        Table with scheduled observations for the given day.
    """
    from astroplan import ObservingBlock
    from astroplan.constraints import TimeConstraint
    from astroplan.scheduling import Schedule
    
    blocks = []
    full_night_start = Time(f"2022-{month}-{day} 19:00") - utcoffset
    full_night_end = Time(f"2022-{month}-{day} 6:00") + 1 * u.day - utcoffset
    full_night = TimeConstraint(full_night_start, full_night_end)
    # Create ObservingBlocks for each filter and target with our time
    # constraint, and durations determined by the exposures needed
    for priority, bandpass in priority_filters:
        # We want each filter to have separate priority (so that target
        # and reference are both scheduled)
        b = ObservingBlock.from_exposures(
            targets,
            priority,
            exp,
            n_exposures,
            read_out,
            configuration={"filter": bandpass},
            constraints=[full_night],
        )
        blocks.append(b)
    noon_before = Time(f"2022-{month}-{day} 12:00") - utcoffset
    noon_after = Time(f"2022-{month}-{day} 12:00") + 1 * u.day - utcoffset
    # Initialize a Schedule object, to contain the new schedule
    day_schedule = Schedule(noon_before, noon_after)
    # Call the schedule with the observing blocks and schedule to schedule the blocks
    scheduler(blocks, day_schedule)
    return day_schedule


def write_plan_to_file(filename, ondate, utcoffset, obs_constraints,
                       target_coord, site_location, site_name):
    """_summary_

    Parameters
    ----------
    filename : _type_
        _description_
    ondate : _type_
        _description_
    utcoffset : _type_
        _description_
    target_coord : astropy.coordinates.SkyCoord
        SkyCoord of the target (must be a single target)
    site_location : _type_
        _description_
    site_name : _type_
        _description_
    """    
    from astroplan import Observer, observability_table
    
    lines = ["UTC, LocalTime, airmass, obs_frac\n"]
    ten_min = 1 * u.hour / 6
    observer = Observer(location=site_location, name=site_name)
    for i in range(23 * 6):
        new_date = ondate + ten_min
        time_range = Time([ondate, new_date])
        obs_table = observability_table(obs_constraints, observer, 
                                        [target_coord], time_range=time_range, 
                                        time_grid_resolution=1*u.minute)
        airmass = target_coord.transform_to(AltAz(obstime=ondate, location=site_location)).secz
        obs_frac = obs_table["fraction of time observable"][0]
        utc = ondate
        lst = ondate + utcoffset
        if obs_table["fraction of time observable"][0] > 0:
            lines.append(f"{utc}, {lst}, {airmass}, {obs_frac}\n")
        ondate = new_date
    with open(filename, "w") as fh:
        fh.writelines(lines)
        
        
def generate_pointing_directions(field_FOV_x, field_FOV_y, detector_FOV_x, detector_FOV_y,
                            finest_level, overlap_frac, center_ra, center_dec, 
                            radius=None):
    """Generate pointing directions given a field of view of the survey field and 
    a detector field of view. The pointing directions are generated in a grid pattern
    to satisfy the dithering pattern. Note that the we can always use a rotation matrix
    to rotate the field if needed, so FOV is really just a 2D angular separation

    Parameters
    ----------
    field_FOV_x : astropy.units.Quantity
        Field of view in the x direction of the survey field.
    field_FOV_y : astropy.units.Quantity
        Field of view in the y direction of the survey field.
    detector_FOV_x : astropy.units.Quantity
        Field of view in the x direction of the detector.
    detector_FOV_y : astropy.units.Quantity
        Field of view in the y direction of the detector.
    finest_level : int
        Finest level of the pointing grid. This is the number of points in the x and y. 
        This is 2^(n-1), where n is the number of levels in the grid.
    overlap_frac : float
        how much of each detector field of view is overlapping. This is the factor
        determines the dithering path.
    center_ra : astropy.units.Quantity
        Center RA of the field.
    center_dec : astropy.units.Quantity
        Center Dec of the field.
    radius : astropy.units.Quantity, optional
        the radius of the circle to constrain the field (has a 5% oversize), by default None

    Returns
    -------
    ra, dec : two 1D-arrays of astropy.units.Quantity
        RA and Dec of the pointing directions.
    """    
    # survey_field_x, survey_field_y = field_FOV_ra, field_FOV_dec
    # detector_FOV_x, detector_y = detector_FOV_ra, detector_FOV_dec
    # overlap for dithering
    overlap_x = overlap_frac * detector_FOV_x
    overlap_y = overlap_frac * detector_FOV_y
    # the finest level the grid is
    # n is essentially the refinement level
    # number of grid is 2^(n-1)
    shift_x = (detector_FOV_x - overlap_x) / 2**(finest_level-1)
    shift_y = (detector_FOV_y - overlap_y) / 2**(finest_level-1)
    num_shift_x = int(np.ceil((field_FOV_x / shift_x).decompose()))
    num_shift_y = int(np.ceil((field_FOV_y / shift_y).decompose()))
    adjusted_field_x = num_shift_x * shift_x
    adjusted_field_y = num_shift_y * shift_y
    # use symmetric padding, i.e. pad the same amount on both sides
    padding_x = adjusted_field_x - field_FOV_x
    padding_y = adjusted_field_y - field_FOV_y
    xi, yi = np.meshgrid(np.linspace(-adjusted_field_x/2, adjusted_field_x/2, num_shift_x), 
                        np.linspace(-adjusted_field_y/2, adjusted_field_y/2, num_shift_y))
    # symmetric padding
    xi -= padding_x/2
    yi -= padding_y/2
    if radius:
        # need to add the circular rejection
        # add a 5% oversize to compensate the loss of image on the side
        in_circle = np.sqrt(xi**2 + yi**2) < radius*1.05
        xi = xi[in_circle]
        yi = yi[in_circle]
    ra = xi.ravel() + center_ra
    dec = yi.ravel() + center_dec
    return ra, dec