from astropy.time import Time
import csv
import datetime
import math
import operator
import requests
import Satellite
import sqlite3
from Site import Site
from Propagation import Observation


def download_tle():
    """
    Reads satellite two-line elements (TLE) from celestrack and saves them as a local text file

    :return: file_string: TLE file name containing the current date in UTC time
    """
    now = datetime.datetime.utcnow()
    file_string = "outputs/" + "TLE " + now.strftime("%m-%d-%Y") + ".txt"

    url = 'https://www.celestrak.com/NORAD/elements/visual.txt'
    r = requests.get(url, allow_redirects=True)
    open(file_string, 'wb').write(r.content)
    return file_string


def read_tle_file(file_name: str, max_sats: int = None) -> list[Satellite]:
    """
    Extracts satellite elements from a TLE file and instantiates a Satellite object for two-line set
    For more information on TLE formats, see https://en.wikipedia.org/wiki/Two-line_element_set

    :param string file_name: Name of the TLE file to read
    :param int max_sats: (Optional: defaults to None) Can be set to reduce the number of satellites read in the case
        of very large TLE files

    :returns: list[Satellite] sat_list: List of Satellite objects
    """
    sat_list = []
    line1_read = False
    line2_read = False
    if max_sats is None:  # If every satellite in the TLE will be read, the limit is set arbitrarily high
        max_sats = 1000000

    with open(file_name, 'r') as tle_file:
        count = 1
        for line in tle_file:
            if count <= max_sats:
                line = ' '.join(line.split())  # Getting rid of duplicated whitespaces found in file
                line = line.split()
                if line:  # Guarding against blank lines in the file
                    if line[0] == '1':  # If this line contains first-line TLE information
                        epoch = line[3]
                        ndot2 = line[4]
                        # ndotdot2 = line[5].replace('-', '.')  # Not needed
                        line1_read = True  # Line 1 has been read

                    elif line[0] == '2':  # If this line contains second-line TLE information
                        incl = line[2]
                        raan = line[3]
                        ecc = '.' + line[4]  # Add a decimal point to the ecc value since TLE files exclude it
                        argp = line[5]
                        mean_anom = line[6]
                        mean_motion = line[7][0:11]  # Slice needed because TLE includes checksum at the end
                        line2_read = True  # Line 2 has been read

                    else:  # If the line is not first-line or second-line information, it's the satellite name
                        sat_name = ' '.join(line)

                    if line1_read and line2_read:
                        # Converting TLE strings into working units
                        ndot2_u, incl_u, raan_u, ecc_u, argp_u, mean_anom_u, mean_motion_u, epoch_jd_u = \
                            tle_data_to_working_units(epoch, ndot2, incl, raan, ecc, argp, mean_anom, mean_motion)

                        # Instantiate new Satellite object
                        new_sat = Satellite.Satellite(
                            sat_name, epoch_jd_u, ndot2_u, incl_u, raan_u, ecc_u, argp_u, mean_anom_u, mean_motion_u
                        )
                        sat_list.append(new_sat)
                        line1_read = False
                        line2_read = False
                        count += 1
    return sat_list


def tle_data_to_working_units(
        epoch: str, ndot2: str, incl: str, raan: str, ecc: str, argp: str, mean_anom: str, mean_motion: str) \
        -> tuple[float, float, float, float, float, float, float, float]:
    """
    Takes raw string params from a TLE text file and converts them into float types with working units

    :param str epoch: Satellite epoch time from TLE (YY + fractional day)
    :param str ndot2: Satellite derivative of mean motion / 2 from TLE (revolutions / day^2)
    :param str incl: Satellite inclination from TLE (deg)
    :param str raan: Satellite right ascension of the ascending node from TLE (deg)
    :param str ecc: Satellite eccentricity from TLE (unitless)
    :param str argp: Satellite argument of perigee from TLE (deg)
    :param str mean_anom: Satellite mean anomaly from TLE (deg)
    :param str mean_motion: Satellite mean motion from TLE (rev / day)

    :returns:
    ndot2: derivative of mean motion/2 (rad/sec^2),
    incl: inclination (rad),
    raan: right ascension of the ascending node (rad),
    ecc: eccentricity (unitless),
    argp: argument of perigee (rad),
    mean_anom: mean anomaly (rad),
    mean_motion: mean_motion (rad/s),
    epoch_jd: satellite epoch (julian day)
    """

    epoch_year = 2000 + int(epoch[0:2])
    epoch_day = float(epoch[2::])
    ndot2 = float(ndot2) * 2 * math.pi / (24 * 60 * 60) ** 2
    ecc = float(ecc)
    incl = float(incl) * math.pi / 180
    raan = float(raan) * math.pi / 180
    argp = float(argp) * math.pi / 180
    mean_anom = float(mean_anom) * math.pi / 180
    mean_motion = float(mean_motion) * 2 * math.pi / (24 * 60 * 60)

    # Epoch day is the fractional day of the year beginning Jan 1 0:00. This is represented as a datetime obj:
    epoch_date = datetime.datetime(year=epoch_year, month=1, day=1) + datetime.timedelta(days=(epoch_day - 1))
    t = Time(epoch_date)  # Creates an astropy time object for conversion to julian days
    epoch_jd = t.jd

    return ndot2, incl, raan, ecc, argp, mean_anom, mean_motion, epoch_jd


def save_observations(site: Site, obs_list: list[Observation], filename: str, new_db=True):
    """
    Saves a list of Observation objects to a sqlite3 database and .csv file.

    The database contains rows formatted for: `satellite name, pass number, jd, local time, range, azimuth, elevation`
    Entries are ordered by the observation start time.

    :param Site site: Site object for the observation location
    :param List[Observation] obs_list: list of all observations to be saved
    :param filename: Name for the output .csv file
    :param new_db: (Optional: defaults to True)
        Setting to True will create a new DB for every new propagation run
        Setting to False will allow the DB to store data for multiple propagations
    """

    db_name = 'outputs/observations.sqlite'
    db = sqlite3.connect(db_name)
    cursor = db.cursor()

    if new_db:
        db.execute("DROP TABLE IF EXISTS "
                   "observations")

    db.execute("CREATE TABLE IF NOT EXISTS "
               "observations ("
               "site TEXT, name TEXT, pass_number REAL, jd REAL, time TEXT, range REAL, azimuth REAL, elevation REAL"
               ")")

    obs_list.sort(key=operator.attrgetter('t0'))  # Sort the observation objects by their start times

    for obs in obs_list:
        datarows = obs.create_observation_entry()

        for row in datarows:
            sat_name, pass_number, jd, sat_range, azimuth, elevation = row

            # Changing the julianday into a readable datetime string
            t = Time(jd, format='jd')
            t.format = 'datetime'
            t += site.utc_offset
            print_time = t.strftime("%m/%d/%Y, %H:%M:%S")

            # Changing units for readability
            azimuth *= 180 / math.pi
            elevation *= 180 / math.pi

            cursor.execute("INSERT INTO observations("
                           "site, name, pass_number, jd, time, range, azimuth, elevation"
                           ") VALUES(?, ?, ?, ?, ?, ?, ?, ?)",
                           (str(site.obs_location), sat_name, pass_number, jd, print_time, sat_range, azimuth,
                            elevation))

    data = cursor.execute("SELECT * FROM observations")

    with open(filename + '.csv', 'w') as obsfile:
        obs_writer = csv.writer(obsfile, delimiter=',')
        obs_writer.writerow([str(site.obs_location)])
        obs_writer.writerow(['Satellite Name', 'Time', 'Range (km)', 'Azimuth (deg)', 'Elevation (deg)'])

        for site, name, pass_number, jd, time, sat_range, azimuth, elevation in data:
            obs_writer.writerow([name, time, round(sat_range, 2), round(azimuth, 2), round(elevation, 2)])

        print("\nObservation details saved to:\n" + filename + '.csv')

    cursor.close()
    db.commit()
    db.close()
