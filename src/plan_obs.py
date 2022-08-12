import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from astropy.units.quantity import Quantity
from astroplan.scheduling import Scheduler


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
    full_night_end = Time(f"2022-{month}-{day} 4:00") + 1 * u.day - utcoffset
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
