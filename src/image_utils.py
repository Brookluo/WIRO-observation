from pathlib import Path


def generate_images_fullpath(input_dir: str, imdir: str, name_format: str, imrange: str):
    """_summary_

    Parameters
    ----------
    input_dir : str
        _description_
    imdir : str
        _description_
    name_format : str
        _description_
    imrange : str
        _description_

    Returns
    -------
    _type_
        _description_
    """        
    # infer the image range from either a contionous range or a list of images
    # should be in the format, 1-10 or 1-10,2-3 or just 1,2,3,5
    full_range = []
    ranges = imrange.split(",")
    for one_range in ranges:
        if "-" in one_range:
            start, end = one_range.split("-")
            full_range.extend(list(range(int(start), int(end)+1)))
        else:
            full_range.append(int(one_range))
    return [Path(input_dir, imdir, name_format.format(i)) for i in full_range]