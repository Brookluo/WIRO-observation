import configparser
from pathlib import Path

from astropy.time import Time

from .image_utils import generate_images_fullpath
from .image_reduction import reduce_images


class WIROObs:
    """_summary_"""

    obs_log_path = None
    images_dict = {}
    done_overscan = False
    done_bias = False
    done_dark = False
    done_flat = False

    def __init__(self, obs_log_path: str):
        self.obs_log_path = obs_log_path
        self.parse_obs_log()
        if self.images_dict["science"]:
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
        if (
            obs_log["bad_images"]["images"] == "NA"
            or obs_log["bad_images"]["images"] == ""
        ):
            self.images_dict["bad_images"] = []
        else:
            # all images myst be their full path for later processing
            self.images_dict["bad_images"] = [
                Path(self.input_dir, obs_log["bad_images"]["imdir"], image.strip())
                for image in obs_log["bad_images"]["images"].split(",")
            ]
        # there are different filters for science and flat correspondingly
        # image_types = ["dark", "bias"]
        sciences = [f"science_{band}" for band in obs_log["science"]["filters"].split(",")]
        flats = [f"flat_{band}" for band in obs_log["flat"]["filters"].split(",")]
        image_types = [*flats, *sciences, "dark", "bias"]
        for imtype in image_types:
            if obs_log[imtype]["imrange"] == "NA" or obs_log[imtype]["imrange"] == "":
                self.images_dict[imtype] = []
            else:
                self.images_dict[imtype] = generate_images_fullpath(
                    self.input_dir,
                    obs_log[imtype]["imdir"],
                    obs_log[imtype]["name_format"],
                    obs_log[imtype]["images"],
                )
    
    def remove_bad_images(self):
        """_summary_
        """        
        # Each file name must be unique for one night of observation!
        # Only science images can be bad!
        self.images_dict["science"] = list(set(self.images_dict["science"]) - set(self.images_dict["bad_images"]))
       
                
    def process_obs_images(self, overscan_sub=True, bias_sub=True, dark_sub=True, flat_div=True, overwrite=False):
        """_summary_"""
        # Keep the hierarchy of the images location after processing?
        bias_img = []
        dark_img = []
        sci_img = []
        if overscan_sub:
           overscan_sub
        if bias_sub:
            self.bias_sub()
        if dark_sub:
            self.dark_sub()
        if flat_div:
            self.flat_div()
