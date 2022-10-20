from pathlib import Path


def assemble_filename(name_format: str, imrange: str):
    """Give a file name format and a range of images, assemble a list of file names.

    Parameters
    ----------
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
    return [name_format.format(i) for i in full_range]


def assemble_fullpath(path, files, stage=None):
    """assemble full paths for a list of files. Users should be able to read
    the full path in the output list.

    Parameters
    ----------
    path : str or Path
        a path to the directory containing the files
    files : list or str
        a list of files or a single file
    stage : str, optional
        which stage images have been processed. Essentially, this is just a suffix
        to the file name, by default None

    Returns
    -------
    list or Path
        a list of full paths to the files or a single Path object
        depends on the input files type
    """    
    if not isinstance(path, Path):
        path_obj = Path(path)
    else:
        path_obj = path
    suffix = ""
    if stage:
        suffix = f"_{stage}"
    if not isinstance(files, list):
        file = Path(files)
        return path_obj / f"{file.stem+suffix}{file.suffix}"
    output_files = []
    for file in files:
        file = Path(file)
        output_files.append(path_obj / f"{file.stem+suffix}{file.suffix}")
    return output_files