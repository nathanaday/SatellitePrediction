import MathCore


class Satellite:
    """
    Class for representing satellite elements.
    Elements include Keplerian elements in addition to rates of change due to J2 and drag effects

    `Attributes`:
        sat_name: Satellite name
        epoch_jd: Satellite's epoch time (from TLE) in julian day format
        ndot2: Second derivative of mean motion / 2 (rad/s^2)
        incl: Orbit inclination (rad)
        raan: Right ascension of the ascending node (rad)
        ecc: Eccentricity (unitless)
        argp: Argument of perigee (rad)
        mean_anom: Mean anomaly (rad)
        mean_motion: Mean motion (rad/s)
        raandot: Right ascension of the ascending node rate of change (rad/s)
        argpdot: Argument of perigee rate of change (rad/s)
        eccdot: Eccentricity rate of change (1/s)

    `Methods`:
        prop(deltat)
            Predicts satellite location given a time of flight, initial orbital elements, and their rates of change
    """

    def __init__(self, sat_name, epoch_jd, ndot2, incl, raan, ecc, argp, mean_anom, mean_motion):
        """
        :param float epoch_jd: Satellite's epoch time (from TLE) in julian day format
        :param float ndot2: Second derivative of mean motion / 2 (rad/s^2)
        :param float incl: Orbit inclination (rad)
        :param float raan: Right ascension of the ascending node (rad)
        :param float ecc: Eccentricity (unitless)
        :param float argp: Argument of perigee (rad)
        :param float mean_anom: Mean anomaly (rad)
        :param float mean_motion: Mean motion (rad/s)
        """

        self.sat_name = sat_name
        self.epoch_jd = epoch_jd
        self.ndot2 = ndot2
        self.incl = incl
        self.raan = raan
        self.ecc = ecc
        self.argp = argp
        self.mean_anom = mean_anom
        self.mean_motion = mean_motion
        self.raandot, self.argpdot, self.eccdot = \
            MathCore.j2_drag(self.incl, self.ecc, self.mean_motion, self.ndot2)  # Rates of change due to J2 and drag

    def prop(self, delta: float) -> tuple[float, float, float, float, float]:
        """
        Predicts satellite location given a time of flight

        :param float delta: Time of flight (seconds)

        :returns:
        Updated mean motion (rad/s),
        Updated eccentricity (unitless),
        Updated right ascension of the ascending node (rad),
        Updated argument of perigee (rad),
        Updated true anomaly (rad)
        """
        mean_motion, ecc, raan, argp, nu = \
            MathCore.coe_update(delta, self.mean_motion, self.ndot2, self.ecc, self.eccdot, self.raan,
                                self.raandot, self.argp, self.argpdot, self.mean_anom)

        return mean_motion, ecc, raan, argp, nu
