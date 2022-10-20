import configparser
from pathlib import Path
import logging
from astropy.time import Time

from image_utils import assemble_filename, assemble_fullpath
from image_reduction import all_overscan_sub_trim, bias_subtract, \
    flat_correct, make_masterflat, make_masterbias


class WIROObs:
    """_summary_"""

    obs_log_path = None
    done_overscan = False
    done_bias = False
    done_dark = False
    done_flat = False
    redo_all = False
    redo_zt = False
    redo_bias = False
    redo_flat = False
    redo_dark = False
    input_dir = None
    output_dir = None
    obs_date = None
    masterbias = None
    masterdark = None


    def __init__(self):
        self.images_dict = {}
        self.masterflats = {}
        self.sci_filters = []
        self.flat_filters = []

    
    def __str__(self) -> str:
        import pprint
        pp = pprint.PrettyPrinter(indent=2, compact=True)
        all_images = f"All images: \n" + pp.pformat(self.images_dict)
        sci_filters = f"Science filters: {self.sci_filters}"
        flat_filters = f"Flat filters: {self.flat_filters}"
        proc = "Processing steps: \n"
        proc += f"Done overscan: {self.done_overscan}\n"
        proc += f"Done bias: {self.done_bias}\n"
        proc += f"Done dark: {self.done_dark}\n"
        proc += f"Done flat: {self.done_flat}"
        dir_info = f"Input directory: {self.input_dir}\n"
        dir_info += f"Output directory: {self.output_dir}"
        return "\n".join([all_images, sci_filters, flat_filters, proc, dir_info])


    def from_lists(self, sci_dict: dict, bias_filelist: list,
                   dark_filelist: list, flat_dict: list, input_dir: str, 
                   output_dir: str, **kwargs):
        """_summary_

        Parameters
        ----------
        sci_filelist : list
            a dictionary of science images, with filter name as key and
            a list of image names as value. Image in the filelist is not
            the full path.
        bias_filelist : list
            _description_
        dark_filelist : list
            _description_
        flat_filelist : list
            _description_
        """
        # this must be {filter: list of files}
        for filt, files in sci_dict.items():
            self.images_dict[f"sci_{filt}"] = files
        for filt, files in flat_dict.items():
            self.images_dict[f"flat_{filt}"] = files
        self.images_dict["bias"] = bias_filelist
        self.images_dict["dark"] = dark_filelist
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.images_dict["bad_images"] = []
        if "obs_date" in kwargs.keys():
            self.obs_date = kwargs["obs_date"]
        if "bad_images" in kwargs.keys():
            self.images_dict["bad_images"] = kwargs["bad_images"]
        self.sci_filters = list(sci_dict.keys())
        self.flat_filters = list(flat_dict.keys())
        self.check_duplicates()
    
    
    def from_obs_log(self, obs_log_path: str):
        """_summary_"""
        self.obs_log_path = obs_log_path
        self.parse_obs_log()
        if self.images_dict["bad_images"]:
            self.remove_bad_images()
        # keep a record of bad images for references
        self._bad_images = self.images_dict["bad_images"]
        del self.images_dict["bad_images"]
        

    def parse_obs_log(self):
        """_summary_"""
        obs_log = configparser.ConfigParser()
        obs_log.optionxform = str
        obs_log.read(self.obs_log_path)
        self.obs_date = Time(obs_log["date"]["obs_date"])
        # this date should be set when the processing step finished
        # self.proc_date = Time(obs_log["date"]["proc_date"])
        self.input_dir = Path(obs_log["locations"]["input_dir"])
        self.output_dir = Path(obs_log["locations"]["output_dir"])
        if (obs_log["bad_images"]["images"] == "NA" 
            or obs_log["bad_images"]["images"] == ""):
            self.images_dict["bad_images"] = []
        else:
            # only the name is necessary, path can be assembled later
            # with input directory information
            self.images_dict["bad_images"] = [image.strip() 
                                              for image in obs_log["bad_images"]["images"].split(",")]
        # there are different filters for science and flat correspondingly
        # image_types = [flat_fileter, sci_filter, "dark", "bias"]
        # Note that if there is science image for a filter, there must be a 
        # flat image for the same filter. However, I assumed a relaxed condition
        self.sci_filters = obs_log["science"]["filters"].split(",")
        self.flat_filters = obs_log["flat"]["filters"].split(",")
        sciences = [f"sci_{band}" for band in self.sci_filters]
        flats = [f"flat_{band}" for band in self.flat_filters]
        image_types = [*flats, *sciences, "dark", "bias"]
        for imtype in image_types:
            if not obs_log.has_section(imtype):
                logging.warning(f"No section for {imtype} in observation log!")
                continue
            if obs_log[imtype]["imrange"] == "NA" or obs_log[imtype]["imrange"] == "":
                self.images_dict[imtype] = []
            else:
                self.images_dict[imtype] = assemble_filename(
                    obs_log[imtype]["name_format"],
                    obs_log[imtype]["imrange"],
                )
        self.check_duplicates()
    
    
    def check_duplicates(self):
        """Check duplicate images in the image list

        Raises
        ------
        RuntimeError
            If there are duplicate images in the image list, 
            raise an Runtime error.
        """        
        from collections import Counter
        all_images = flatten(list(self.images_dict.values()))
        counts = Counter(all_images)
        dup = False
        for image, count in counts.items():
            if count == 1:
                continue
            print(f"Duplicate image: {image}\n")
            print("Occurences:\n")
            for key, imlist in self.images_dict.items():
                if image in imlist:
                    print(f"{key}\n")
            dup = True
        if dup:
            raise RuntimeError("Duplicate images found!")


    def remove_bad_images(self):
        """_summary_
        """        
        # Each file name must be unique for one night of observation!
        # Only science images can be bad!
        for filt in self.sci_filters:
            self.images_dict[f"sci_{filt}"] = [
                image for image in self.images_dict[f"sci_{filt}"]
                if image not in self.images_dict["bad_images"]
            ]


    def overscan_sub_trim(self):
        """_summary_"""
        all_images = flatten(list(self.images_dict.values()))
        files_fullpath = assemble_fullpath(self.input_dir, all_images)
        all_overscan_sub_trim(files_fullpath, self.output_dir, polyfit='cheb', order=3,
                              overwrite=self.redo_all|self.redo_zt)
        self.done_overscan = True
        

    def bias_sub(self):
        """_summary_"""
        if not self.images_dict["bias"]:
            # logging.error("No bias images for bias subtraction!")
            raise ValueError("No bias images for bias subtraction!")
        if self.done_overscan:
            last_stage = "z"
            input_dir = self.output_dir
        else:
            last_stage = ""
            input_dir = self.input_dir
        bias_img = assemble_fullpath(input_dir, self.images_dict["bias"], last_stage)
        all_images = flatten(list(self.images_dict.values()))
        # remove bias images from all images
        all_images = list(set(all_images) - set(self.images_dict["bias"]))
        all_fullpath = assemble_fullpath(input_dir, all_images, last_stage)
        make_masterbias(bias_img, self.output_dir, overwrite=self.redo_all|self.redo_bias)
        self.masterbias = Path(self.output_dir, "masterbias.fits")
        bias_subtract(all_fullpath, self.masterbias, self.output_dir, 
                        overwrite=self.redo_all|self.redo_bias)
        self.done_bias = True
 

    def dark_sub(self):
        pass
    

    def flat_div(self):
        """_summary_
        """        
        if not self.done_overscan:
            input_dir = self.input_dir
            last_stage = ""
        elif not self.done_bias:
            input_dir = self.output_dir
            last_stage = "z"
        else:
            input_dir = self.output_dir
            last_stage = "zb"
        path_masterflat_tmpl = "masterflat_{0}_clip_med_weighted_count.fits"
        # path_masterflat_tmpl = "masterflat_{0}_norm.fits"
        for band in self.flat_filters:
            make_masterflat(
                assemble_fullpath(input_dir, self.images_dict[f"flat_{band}"], last_stage), 
                self.output_dir, band, overwrite=self.redo_all|self.redo_flat,
                save_plots=True
                )
            self.masterflats[band] = Path(self.output_dir, path_masterflat_tmpl.format(band))
        all_science = [v for k, v in self.images_dict.items() if k.startswith("sci_")]
        all_sci_path = assemble_fullpath(input_dir, flatten(all_science), last_stage)
        flat_correct(all_sci_path, self.masterflats, self.output_dir, overwrite=self.redo_all|self.redo_flat)
        self.done_flat = True


    def process_obs_images(self, overscan_sub=True, bias_sub=True, dark_sub=True, flat_div=True):
        """_summary_"""
        if overscan_sub:
            logging.info("Doing overscan subtraction and trimming...")
            self.overscan_sub_trim()
        if bias_sub:
            logging.info("Doing bias subtraction...")
            self.bias_sub()
        if dark_sub:
            logging.info("Doing dark subtraction...")
            self.dark_sub()
        if flat_div:
            logging.info("Doing flat correction...")
            self.flat_div()
        logging.info("Image reduction finished!")


def flatten(l):
    return [item for sublist in l for item in sublist]