import config
import datetime
import math
import numpy as np
import pytz
import requests
from timezonefinder import TimezoneFinder


class Site:
    """
    Class for representing observations sites

    `Attributes`:
        obs_location: geopy geocode location object
        name: (optional, defaults to None) Site name
        lon: Site longitude (deg)
        lat: Site latitude (deg)
        timezone: Site timezone
        utc_offset: Datetime timedelta for the difference between local time and UTC time at site location
        altitude: Site altitude (m)

    `Methods`:
        set_name(name)
            Setter for assigning a name to the Site

        find_info_from_location()
            Extracts the longitude and latitude from the geopy geolocation and assigns them to the class attributes
            Finds the timezone at the site location based on longitude and latitude with TimeZoneFinder()

        find_elevation()
            Uses the open-elevation API to make a query with the site's lon and lat, then sets the result as the
            site's altitude attribute.

        find_utc_offset()
            Finds the difference between local time and utc time as a datetime timedelta, then sets the timedelta
            as a site attribute

        site_vector(site_lst)
            Calculates and returns the ECI position vector of the site
    """

    def __init__(self, obs_location):
        """
        :param obs_location: geopy geocode for site location
        """
        self.obs_location = obs_location
        self.name = None
        self.lon = None
        self.lat = None
        self.timezone = None
        self.utc_offset = None
        self.altitude = None
        self.find_info_from_location()
        self.find_elevation()
        self.find_utc_offset()

    def set_name(self, name: str):
        """
        Setter for assigning a name to the Site
        :param str name: Name for the site
        """
        self.name = name

    def find_info_from_location(self):
        """
        Extracts the longitude and latitude from the geopy geolocation and assigns them to the class attributes
        Finds the timezone at the site location based on longitude and latitude with TimeZoneFinder()
        """
        self.lon = self.obs_location.longitude
        self.lat = self.obs_location.latitude
        self.timezone = TimezoneFinder().timezone_at(lng=self.lon, lat=self.lat)

    def find_elevation(self):
        """
        Uses the open-elevation API to make a query with the site's lon and lat, then sets the result as the
        site's altitude attribute.
        """
        try:
            query = ('https://elevation-api.io/api/elevation?points=({},{})&key={}'.format(
                    self.lat, self.lon, config.elevation_api))
            r = requests.get(query).json()
            self.altitude = r["elevations"][0]["elevation"]
        except KeyError:
            print('-'*10 + " WARNING " + '-'*10)
            print("There was a problem retrieving elevation from elevation-api.\n"
                  "It's possible the API hasn't been added, or was added incorrectly.\n"
                  "See README for instructions on how to easily link the API with config.py\n"
                  "The program can continue with a default site altitude of 0km, but note this may produce skewed "
                  "results in locations at significant altitude")
            print('-'*10 + " WARNING " + '-'*10)
            self.altitude = 0

    def find_utc_offset(self):
        """
        Finds the difference between local time and utc time as a datetime timedelta, then sets the timedelta
        as a site attribute
        """
        site_now = datetime.datetime.now(pytz.timezone(self.timezone))
        delta = site_now.utcoffset().total_seconds() / 60 / 60
        self.utc_offset = datetime.timedelta(hours=delta)

    def site_vector(self, site_lst: float) -> np.ndarray:
        """
        Calculates and returns the ECI position vector of the site

        :param float site_lst: local sidereal time (rad)
        :return: numpy ndarray representing the ECI position vector for the site (km)
        """

        # Global values:
        earth_radius = 6378.137
        earth_eccentricity = 0.0818191908

        # Convert to working units
        lat = self.lat * math.pi / 180
        elevation = self.altitude / 1000

        x = abs(earth_radius / math.sqrt(1 - earth_eccentricity ** 2 * (math.sin(lat)) ** 2)
                + elevation) * math.cos(lat)

        z = abs(earth_radius * (1 - earth_eccentricity ** 2) /
                math.sqrt(1 - earth_eccentricity ** 2 * (math.sin(lat)) ** 2) + elevation) * math.sin(lat)

        site_vector = np.array([[x * math.cos(site_lst)], [x * math.sin(site_lst)], [z]])  # Column vector

        return site_vector
