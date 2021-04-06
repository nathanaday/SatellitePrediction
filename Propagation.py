import Satellite
from Site import Site
import MathCore
import math
from datetime import datetime
from astropy.time import Time


class Observation:
    """
    Class for representing a satellite observation opportunity.
    Because propagating observations is done iteratively through time, the class stores data in lists that are appended
    for each iteration of the observation opportunity

    `Attributes`:
        sat_name: Satellite name
        pass_number: Indicates a place order among other observations; serves as a unique identifier
        time_jd: List of times in julian day format for each iteration of this observation
        sat_range: List of satellite ranges for each iteration of this observation (km)
        azimuth:  List of satellite azimuth angles for each iteration of this observation (rad)
        elevation: List of satellite elevation angles for each iteration of this observation (rad)
        t0: The start time of this observation (julian day)

    `Methods`:
        add_data(new_timejd, new_satrange, new_azimuth, new_elevation)
            Appends the class's data lists (time, range, azimuth, and elevation) with new data from a pass iteration

        create_observation_entry()
            Consolidates the observation data lists into a list of tuples formatted for data output
    """

    def __init__(self, sat_name, pass_number):
        """
        :param str sat_name: Satellite name
        :param int pass_number: Indicates a place order among other observations; serves as a unique identifier
        """
        self.sat_name = sat_name
        self.pass_number = pass_number
        self.time_jd = []
        self.sat_range = []
        self.azimuth = []
        self.elevation = []
        self.t0 = None

    def add_data(self, new_timejd: float, new_satrange: float, new_azimuth: float, new_elevation: float):
        """
        Appends the class's data lists (time, range, azimuth, and elevation) with new data from a pass iteration

        :param float new_timejd: Time at iteration in julian day format
        :param float new_satrange: Sat range at iteration (km)
        :param float new_azimuth: Sat azimuth angle at iteration (rad)
        :param float new_elevation: Sat elevation at iteration (rad)
        """
        if not [x for x in (new_timejd, new_satrange, new_azimuth, new_elevation) if x is None]:  # If none are empty
            self.time_jd.append(new_timejd)
            self.sat_range.append(new_satrange)
            self.azimuth.append(new_azimuth)
            self.elevation.append(new_elevation)

            self.t0 = min(self.time_jd)  # Updating the start time

    def create_observation_entry(self) -> list[tuple[str, int, float, float, float, float]]:
        """
        Consolidates the observation data lists into a list of tuples formatted for data output

        :return:
        datarows: A list of tuples with the format:
        (satellite name, pass number, time jd, range km, azimuth rad, elevation rad)
        """
        datarows = []
        if len({len(i) for i in
                [self.time_jd, self.sat_range, self.azimuth, self.elevation]}) == 1:  # If all lists are same length
            for i, time_jd in enumerate(self.time_jd):
                row = (self.sat_name, self.pass_number, time_jd, self.sat_range[i], self.azimuth[i], self.elevation[i])
                datarows.append(row)
            return datarows
        else:
            print("Not all the data rows have the same length")


def propagate_satellites(sats: list[Satellite.Satellite], start_time: datetime, stop_time: datetime, site: Site) \
        -> list[Observation]:
    """
    Propagates satellite orbits and returns instances where the satellite is visible from the observation location

    :param list[Satellite] sats: A list of satellite objects to be propagated (min size 1)
    :param datetime start_time: Datetime object for the beginning of the desired observation window
    :param datetime stop_time: Datetime object for the end of the desired observation window
    :param Site site: Reference to the Site object to indicate the location of observation

    :return: obs: A list of Observation objects for each instance of a visible pass
    """

    obs_list = []
    pass_counter = 0
    obs_last_pass = False

    lon = site.lon * math.pi / 180  # Convert to working units

    # Convert the  starting and stopping datetimes to Astropy Time objects with julian day format
    t_start = Time(start_time.strftime("%Y-%m-%d %H:%M:%S"), format='iso', scale='utc')
    t_start.format = 'jd'
    t_stop = Time(stop_time.strftime("%Y-%m-%d %H:%M:%S"), format='iso', scale='utc')
    t_stop.format = 'jd'

    # julian day start and stop values as python floats
    t_stop_f = t_stop.to_value('jd').item()
    t_start_f = t_start.to_value('jd').item()

    for i, sat in enumerate(sats):
        print("({0} / {1}) Propagating orbit for {2}".format(i+1, len(sats), sat.sat_name))
        obs = None  # Initializing observation variable

        # Seconds between sat epoch and start time:
        delta_t1s = (t_start_f - sat.epoch_jd) * 24 * 60 * 60

        working_jd_f = t_start_f
        while working_jd_f <= t_stop_f:

            new_mean_motion, new_ecc, new_raan, new_argp, new_nu = sat.prop(delta_t1s)  # Propagate orbital elements

            # Greenwich sidereal time and local sidereal time
            gst = MathCore.find_gst(working_jd_f)
            lst = gst + lon

            r_site_vec = site.site_vector(lst)  # Site position vector (ECEF frame)
            r_sat_vec = MathCore.sat_vec(
                new_mean_motion, new_ecc, new_nu, new_argp, sat.incl, new_raan)  # Sat position vector (ECEF frame)

            sat_range, azimuth, elevation = \
                MathCore.range_elevation_azimuth(r_sat_vec, r_site_vec, site, lst)  # ECEF to celestial coord frame

            # Checking if the satellite is visible during this iteration...
            if MathCore.is_visible(working_jd_f, r_sat_vec, r_site_vec, elevation):
                obs_this_pass = True
                if not obs_last_pass:  # If the last iteration wasn't visible, i.e. first iteration of new obs
                    pass_counter += 1
                    obs = Observation(sat.sat_name, pass_counter)  # Instantiate new Observation object
                    obs_list.append(obs)
                obs.add_data(working_jd_f, sat_range, azimuth, elevation)
            else:
                obs_this_pass = False

            obs_last_pass = obs_this_pass  # Setting `this pass` as `last pass` for next iteration
            working_jd_f += 100 / 60 / 60 / 24  # Adds 100 seconds in `days`
            delta_t1s = (working_jd_f - sat.epoch_jd) * 24 * 60 * 60  # New gap between working time and sat epoch

    print("\nThere are {} viewing opportunities tonight".format(pass_counter))
    return obs_list
