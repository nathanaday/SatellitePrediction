import math
import numpy as np
from fractions import Fraction
from Site import Site


def is_numpy(a):
    """
    Raises an exception if an expected numpy array is not a numpy array

    :param a: Value to check
    :return: True if `a` is a numpy array
    """
    if isinstance(a, np.ndarray):
        return True
    else:
        raise Exception("Expected a numpy array for vector input")


def vecangle(vector1: np.ndarray, vector2: np.ndarray) -> float:
    """
    Calculates the angle in radians between two vectors

    :param np.ndarray vector1: The first vector
    :param np.ndarray vector2: The second vector
    :returns: The angle between the first and second vectors in radians
    """
    is_numpy(vector1)
    is_numpy(vector2)

    unit_v1 = vector1 / np.linalg.norm(vector1)
    unit_v1 = unit_v1.flatten()
    unit_v2 = vector2 / np.linalg.norm(vector2)
    unit_v2 = unit_v2.flatten()

    angle = np.arccos(np.dot(unit_v1, unit_v2))
    return angle


def rotate_axis(matrix: np.ndarray, axis: int, angle: float) -> np.ndarray:
    """
    Allows a 3x1 matrix to be rotated about one axis by an angle amount.

    :param np.ndarray matrix: (3x1) The input matrix before rotation
    :param int axis: The axis of rotation (1, 2, or 3, corresponding to the dimensions of the matrix)
    :param float angle: The angle of rotation in radians

    :returns: matrix: (3x1) Numpy matrix after rotation
    """
    is_numpy(matrix)
    el1, el2, el3 = None, None, None  # Initializing elements

    try:
        if axis == 1:
            el1 = matrix[0][0]
            el2 = matrix[1][0] * math.cos(angle) + matrix[2][0] * math.sin(angle)
            el3 = -matrix[1][0] * math.sin(angle) + matrix[2][0] * math.cos(angle)
        if axis == 2:
            el1 = matrix[0][0] * math.cos(angle) - matrix[2][0] * math.sin(angle)
            el2 = matrix[1][0]
            el3 = matrix[0][0] * math.sin(angle) + matrix[2][0] * math.cos(angle)
        if axis == 3:
            el1 = matrix[0][0] * math.cos(angle) + matrix[1][0] * math.sin(angle)
            el2 = -matrix[0][0] * math.sin(angle) + matrix[1][0] * math.cos(angle)
            el3 = matrix[2][0]

        matrix = np.array([[el1], [el2], [el3]])

    except (ValueError, TypeError):
        print("Invalid matrix, axis, or angle in rotate_axis()")
    return matrix


def newtons_method(mean_anomaly: float, eccentricity: float) -> float:
    """
    Uses newton's method (also know as Newton-Raphson method) to iteratively find better approximations
    of an orbit's eccentric anomaly.

    :param float mean_anomaly: The orbit's mean anomaly (radians)
    :param float eccentricity: The orbit's eccentricity (unitless)

    :returns: eccentric_anom: best approximation for eccentric anomaly (radians)
    """

    eccentric_anom = mean_anomaly  # Arbitrary first guess for eccentric-anomaly
    difference = 10  # Arbitrary difference to begin with

    while difference > 10 ** -10:
        e_new = eccentric_anom + (mean_anomaly - eccentric_anom - eccentricity * math.sin(eccentric_anom)) / (
                1 - eccentricity * math.cos(eccentric_anom))
        difference = abs(e_new - eccentric_anom)
        eccentric_anom = e_new
    return eccentric_anom


def truncate_jd(jd: float) -> float:
    """
    Takes only the decimal value of any jd and adds 1 if the UTC time is past noon, allowing JD times to be compared
    as fractional days

    :param float jd: Julian day float
    :returns: fractional day float
    """
    decimal_jd = jd % 1
    if decimal_jd < 0.5:
        decimal_jd += 1
    return decimal_jd


def j2_drag(incl: float, ecc: float, mean_motion: float, ndot2: float) -> tuple[float, float, float]:
    """
    Calculates rates of chance for right ascension of the ascending node, argument of perigee, and
    eccentricity due to the J2 effect and aerodynamic drag

    :param float incl: Initial orbit inclination (radians)
    :param float ecc: Initial orbit eccentricity (unitless)
    :param float mean_motion: Initial mean-motion (rad/sec)
    :param float ndot2: Initial mean-motion rate divided by 2 (rad/sec^2)

    :returns:
    raan rate of change (rad/sec),
    argp rate of change (rad/sec),
    eccentricity rate of change (1/sec)
    """

    # Define constants
    j2 = 1.08263 * 10 ** -3
    r0 = 6378.137

    # Orbital parameter p0:
    a = (398600.5 / (mean_motion ** 2)) ** Fraction(1 / 3)
    p0 = a * (1 - ecc ** 2)

    # New mean motion with J2
    n = mean_motion * (1 + 3 / 2 * j2 * (r0 / p0) ** 2 * math.sqrt(1 - ecc ** 2) * (1 - 3 / 2 * math.sin(incl) ** 2))

    # COE rates of change
    raandot = (-3 / 2 * j2 * (r0 / p0) ** 2 * math.cos(incl)) * n
    argpdot = (3 / 2 * j2 * (r0 / p0) ** 2 * (2 - 5 / 2 * (math.sin(incl)) ** 2)) * n
    eccdot = -(2 * (1 - ecc) * ndot2 * 2) / (3 * mean_motion)

    return raandot, argpdot, eccdot


def coe_update(delta_t: float, mean_motion: float, ndot2: float, ecc: float, eccdot: float, raan: float, raandot: float,
               argp: float, argpdot: float, mean_anomaly: float) -> tuple[float, float, float, float, float]:
    """
    Updates an orbit's central orbital elements (COEs) through a specified time of flight given their rates of change

    :param float delta_t: The amount of time since the initial orbital elements and the updated set (seconds)
    :param float mean_motion: Orbital mean motion at t0 (rad/sec)
    :param float ndot2: Orbital rate of change of mean motion divided by 2 at t0 (rad/sec^2)
    :param float ecc: Eccentricity at t0 (unitless)
    :param float eccdot: Eccentricity rate of change due to perturbations (1/sec)
    :param float raan: Right ascension of the ascending node at t0 (rad)
    :param float raandot: Right ascension of the ascending node rate of change due to perturbations (rad/sec)
    :param float argp: Argument of perigee at t0 (rad)
    :param float argpdot: Argument of argument of perigee rate of change (rad/sec)
    :param float mean_anomaly: Orbital mean anomaly parameter (rad)

    :returns:
    Updated mean motion (rad/s),
    Updated eccentricity (unitless),
    Updated right ascension of the ascending node (rad),
    Updated argument of perigee (rad),
    Updated true anomaly (rad)
    """

    # Future central orbital elements
    new_mean_motion = mean_motion + ndot2 * 2 * delta_t
    new_ecc = ecc + eccdot * delta_t
    new_raan = raan + raandot * delta_t
    new_argp = argp + argpdot * delta_t

    # Future mean anomaly
    future_mean_anom = mean_anomaly + mean_motion * delta_t + ndot2 * delta_t ** 2

    # Future eccentric anomaly
    eccentric_anom = newtons_method(future_mean_anom, new_ecc)

    # The new true anomaly is the inverse cosine of two parameters, p1 and p2:
    p1 = (math.cos(eccentric_anom) - new_ecc)
    p2 = (1 - new_ecc * math.cos(eccentric_anom))
    p12 = p1 / p2
    new_nu = math.acos(p12)

    # Quadrant check on true anomaly
    if eccentric_anom > math.pi:
        new_nu = 2 * math.pi - new_nu

    return new_mean_motion, new_ecc, new_raan, new_argp, new_nu


def range_elevation_azimuth(sat_vec: np.ndarray, r_site_vector: np.ndarray, site_ref: Site, lst: float):
    """
    Finds celestial coordinates (azimuth, elevation) of a satellite as seen from an observation site, as well as range

    :param np.ndarray sat_vec: Position vector of the satellite in the ECI system (km)
    :param np.ndarray r_site_vector: Position vector of the observation site in ECI (km)
    :param Site site_ref: Reference to the observation site object
    :param float lst: Local sidereal time (rad)

    :returns:
    rho: Satellite range (km),
    azimuth: Satellite azimuth angle (rad)
    elevation: Satellite elevation angle (rad)

    """
    is_numpy(sat_vec)

    lat = site_ref.lat * math.pi / 180  # Convert to working units

    # Subtracting satellite position vector from site vector gives the range between the site and satellite
    rho_ijk = sat_vec - r_site_vector
    rho_ijk_shaped = rho_ijk.reshape(-1, 1)

    # The scalar range is just the magnitude of the range vector
    rho = np.linalg.norm(rho_ijk_shaped)

    # Defining parameters b, y to make matrix operation more readable
    b = math.pi / 2 - lat
    y = lst
    transform1 = np.array([[math.cos(b), 0, -math.sin(b)], [0, 1, 0], [math.sin(b), 0, math.cos(b)]])
    transform2 = np.array([[math.cos(y), math.sin(y), 0], [-math.sin(y), math.cos(y), 0], [0, 0, 1]])

    # Finding rho in the south-east-zenith (SEZ) coordinate frame
    mat1 = np.matmul(transform1, transform2)
    rho_sez = np.matmul(mat1, rho_ijk)

    # Finding the elevation using rho_sez and the range magnitude
    elevation = math.asin(rho_sez[2][0] / rho)

    if elevation > math.pi:  # If elevation is greater than 180 degs (over the horizon)
        elevation = 2 * math.pi - elevation

    # Finding azimuth using rho_sez, the range magnitude, and elevation
    azimuth = math.pi - math.acos(rho_sez[0][0] / (rho * math.cos(elevation)))

    # The above calculation only produces an azimuth between 0 and pi rads. We can do a quadrant check to compensate:
    if rho_sez[1][0] < 0:  # If the east component of the range vector < 0 (azimuth between pi and 2pi rads)
        azimuth = 2 * math.pi - azimuth  # Azimuth moves to the next quadrant

    return rho, azimuth, elevation


def sat_vec(mean_motion: float, ecc: float, nu: float, argp: float, inc: float, raan: float) -> np.ndarray:
    """
    Finds the ECI position vector of a satellite given its orbital elements

    :param mean_motion: orbit's mean motion (rad/s)
    :param ecc: orbit's eccentricity (unitless)
    :param nu: orbit's true anomaly (rad)
    :param argp: orbit's argument of perigee (rad)
    :param inc: orbit's inclination (rad)
    :param raan: orbit's right ascension of the ascending node (rad)

    :return: r_vec_ECI: satellite position vector in the ECI frame (km)
    """
    expo = Fraction(1 / 3)
    semimajoraxis = (398600.5 / mean_motion ** 2) ** expo

    # Magnitude of the satellite's range vector using the two-body equation of motion
    r_mag = semimajoraxis * (1 - ecc ** 2) / (1 + ecc * math.cos(nu))

    # The range vector in the perifocal reference frame
    r_vec_pqw = np.array([[r_mag * math.cos(nu)], [r_mag * math.sin(nu)], [0]])

    # Three matrix rotations are required to move from PQW frame to the ECI frame
    rot1 = rotate_axis(r_vec_pqw, 3, -argp)  # Rotate about axis 3 by (-) argument of perigee
    rot2 = rotate_axis(rot1, 1, -inc)  # Rotate about axis 1 by (-) inclination
    r_vec_ECI = rotate_axis(rot2, 3, -raan)  # Rotate about axis 3 by (-) right asc. of asc. node

    return r_vec_ECI


def sun_vector(julian_day: float) -> tuple[float, np.ndarray]:
    """
    Finds the ECI position vector of the sun on a given julian day.
    Note: this is an approx. but has a precision of <0.1 deg, and is much faster than high-precision alternatives

    Equations from USNO approximate solar coordinates:
    https://web.archive.org/web/20181115153648/http://aa.usno.navy.mil/faq/docs/SunApprox.php
    Also detailed in this wiki:
    https://en.wikipedia.org/wiki/Position_of_the_Sun#Ecliptic_coordinates

    :param float julian_day: julian day

    :return:
    sun_range: distance from earth to sun (km)
    sun_vec: ECI position vector of the sun (km)
    """

    # Constant
    obliquity = 23.45  # Deg

    days_since_2000 = julian_day - 2451545.0
    mean_lon = 280.460 + 0.9856474 * days_since_2000  # Sun mean longitude (deg)
    mean_lon = math.fmod(mean_lon, 360)
    if mean_lon < 0:
        mean_lon += 360

    mean_anomaly = 357.528 + 0.9856003 * days_since_2000  # Sun mean anomaly (deg)
    mean_anomaly = math.fmod(mean_anomaly, 360)
    if mean_anomaly < 0:
        mean_anomaly += 360

    ecliptic_lon = mean_lon + 1.915 * math.sin(mean_anomaly * math.pi / 180) + \
                   0.02 * math.sin(2 * mean_anomaly * math.pi / 180)  # Sun ecliptic longitude (deg)
    if ecliptic_lon < 0:
        ecliptic_lon += 360

    mean_lon *= math.pi / 180  # deg -> radians
    mean_anomaly *= math.pi / 180  # deg -> radians
    ecliptic_lon *= math.pi / 180  # deg -> radians
    obliquity *= math.pi / 180  # deg -> radians

    sun_range = (1.00014 - 0.01671 * math.cos(mean_anomaly) - 0.00014 * math.cos(2 * mean_anomaly)) * (1.496 * 10 ** 8)

    sun_vec = np.array([
        [sun_range * math.cos(ecliptic_lon)],
        [sun_range * math.cos(obliquity) * math.sin(ecliptic_lon)],
        [sun_range * math.sin(obliquity) * math.sin(ecliptic_lon)]
    ])

    # return beta, ecliptic_lon
    return sun_range, sun_vec


def is_visible(jd: float, sat_vec: np.ndarray, site_vec: np.ndarray, elevation: float) -> bool:
    """
    Determines if a satellite is visible given its position and an observation site

    :param jd: Julian day at time of visibility check
    :param sat_vec: ECI position vector of the satellite (km)
    :param site_vec: ECI position vector of the site (km)
    :param elevation: Satellite's elevation as seen from observation site (rad)

    :returns: boolean corresponding to satellite visibility
    """
    is_numpy(sat_vec)

    is_lit = False
    is_dark = False
    in_view = False

    # Checking if the satellite is illuminated by the sun:
    sun_range, sun_vec = sun_vector(jd)
    beta = vecangle(sat_vec, sun_vec)  # Angle between satellite vector and sun vector
    x = np.linalg.norm(sat_vec) * math.sin(math.pi - beta)  # Perpendicular component of the satellite wrt sun
    if abs(x) > 6378.137:  # If the satellite's x-component is greater than earth's radius (e.g. outside the shadow)
        is_lit = True

    # Checking if the site is dark
    alpha = vecangle(site_vec, sun_vec)  # Angle between obs location and sun (0 ~= noon, 180 ~= midnight, etc.)
    if 102 <= alpha * (180 / math.pi) <= 258:  # Defines "dark" as about an hour after sunset/ an hour before sunrise
        is_dark = True

    # Checking if the satellite is visible above the horizon (with 6 degrees to compensate for city/landscape)
    if (elevation * 180 / math.pi) > 6:
        in_view = True

    return is_lit and in_view and is_dark


def find_gst(jd: float) -> float:
    """
    Finds the greenwich sidereal time given a julian day

    :param float jd: Julian day

    :return: greenwich sidereal time (rad)
    """
    # Finding the julian day at UTC midnight
    jd_min = int(jd) - 0.5
    jd_max = int(jd) + 0.5
    if jd >= jd_max:
        jd_fix = jd_max
    else:
        jd_fix = jd_min

    #  Standard equation for Greenwich sidereal time at midnight (sidereal seconds)
    t_u = (jd_fix - 2451545) / 36525
    h0 = 24110.54841 + 8640184.812866 * t_u + 0.093104 * (t_u ** 2) - (6.2 * 10 ** -6) * (t_u ** 3)
    sidereal_rate = 1.00273790935 + 5.9 * (10 ** -11) * t_u

    # The number of seconds since UTC midnight
    utc_seconds = (jd - jd_fix) * 86400

    # Greenwich sidereal time (radians)
    gst = (h0 + sidereal_rate * utc_seconds) * (2 * math.pi / 86400)
    gst = gst % (2 * math.pi)
    return gst
