import config
import Data
import datetime
from geopy.geocoders import Nominatim
import os
from pathlib import Path
from Propagation import propagate_satellites
from Site import Site
import sys


def main():
    """
    Reads satellite elements from TLE, takes user input, then propagates all satellites and saves all observations.

    `Steps`

    1. Downloads two-line element (TLE) file and read elements into a list of Satellite objects
    2. Sets observation time frame
        By default, time begins when the program is run and ends 24 hours later
    3. Instantiates observation site.
        The user inputs a general location and geopy matches an exact geocode location
        The user has the option to name the site so that it appears in the .csv file output, but this optional
    4. Propagates satellites
        The propagator predicts the orbit of each satellite in the Sat list and finds opportunities where the satellite
        is visible from the observation location after sunset and before sunrise.
    5. Saves observations
        All predicted opportunities are saved to a sqlite3 database and a .csv file
    """

    # Confirming the `outputs` directory exists
    cwd_path = os.getcwd()
    Path(cwd_path + '/outputs').mkdir(parents=False, exist_ok=True)

    file_string = Data.download_tle()  # Downloads satellite elements from celestrack, gets file name in return
    new_sat_list = Data.read_tle_file(file_string)  # Reads orbital elements, returns list of satellite objects

    start_time = datetime.datetime.utcnow()
    stop_time = start_time + datetime.timedelta(days=1)

    site_name = input("Enter observation location (city, street, or nearby landmark...):\n")
    try:
        geolocator = Nominatim(user_agent="SatellitePrediction")
        location = geolocator.geocode(site_name)
        new_site = Site(location)  # Instantiating Site object

    except AttributeError:
        print("{} didn't match to a longitude or latitude. "
              "Try entering the location again with a different representation\n "
              "(Tip: check spelling, or try a zip code.)".format(site_name))
        sys.exit(1)

    print()
    print(new_site.obs_location)
    print("Lon: {}\nLat: {}\nElv: {}m\nTZ: {}\n".format(
        new_site.lon, new_site.lat, new_site.altitude, new_site.timezone))

    new_site.set_name(input("(Optional: Enter to skip) "
                            "Name this site for the output file\tEx: Home, Jakarta, Pikes Peak\n"))

    if new_site.name != '':  # If there is a full site name, it needs to be cleaned for the filename
        new_site.name = new_site.name.rstrip().replace(',', '_').replace('.', '_').replace(' ', '_')
        new_site.name = '(' + new_site.name + ')'
    local_date = (start_time + new_site.utc_offset).strftime("%m-%d-%Y")
    obs_file = "outputs/OBS {}{}".format(local_date, new_site.name)

    # Satellite propagation:
    obs_list = propagate_satellites(new_sat_list, start_time, stop_time, new_site)
    Data.save_observations(new_site, obs_list, obs_file)


if __name__ == "__main__":
    main()
