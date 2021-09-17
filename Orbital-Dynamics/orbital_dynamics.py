import numpy as np


def periapsis_to_true_anomaly(eccentricity, angular_momentum, true_anomaly, mu=398600):
    """
    :param eccentricity: The eccentricity of an elliptical orbit
    :param angular_momentum: The angular momentum of an elliptical orbit
    :param true_anomaly: The true anomaly of the desired point
    :param mu: The gravitational parameter of the body, set to Earth by default
    :return: Returns the time from the periapsis to the inputted true anomaly
    """

    orbit_period = (2 * np.pi / mu ** 2) * (angular_momentum / np.sqrt(1 - eccentricity ** 2)) ** 3

    eccentric_anomaly = 2 * np.arctan(np.sqrt((1 - eccentricity) / (1 + eccentricity) * np.tan(true_anomaly / 2)))
    mean_anomaly = eccentric_anomaly - eccentricity * np.sin(eccentric_anomaly)

    time = (mean_anomaly / (2 * np.pi)) * orbit_period

    return time
