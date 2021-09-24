import numpy as np

def gauss_variational_impulse(orbital_elements, mu, true_anomaly, delta_v):

    """
    The impulse form of Gauss' variational equations

    Args:
        orbital_elements: List of orbital elements, in the format
            [semi_major_axis, eccentricity, inclination, right_ascension, argument_perigee, true_anomaly]
            Angular units must be in radians

        mu: The gravitational parameter of the body, defined as mu = gravitational constant * body mass

        true_anomaly: The true anomaly of the orbiting body, in radians

        delta_v: List of the changes in velocity, in the format [delta_v_radial, delta_v_normal, delta_v_angular]
            delta_v_angular must be in rad/s

    Returns:
        A list containing the changes in orbital parameters in the form:
            [delta_a, delta_e, delta_i, delta_asc, delta_arg]
            Angular units are in radians
    """

    a, e, i, asc, arg = orbital_elements
    delta_v_radial, delta_v_normal, delta_v_angular = delta_v

    sqrt_amu = np.sqrt(a*(1 - e**2)/mu)

    delta_a = (2*a**2 / np.sqrt(mu*a*(1-e**2))) * (e*np.sin(true_anomaly*delta_v_radial) + (1 + e*np.cos(true_anomaly))*delta_v_angular)

    delta_e = sqrt_amu * (np.sin(true_anomaly*delta_v_radial) + delta_v_angular(2*np.cos(true_anomaly) + e*(1 + np.cos(true_anomaly)**2))/())

    delta_i = sqrt_amu * delta_v_normal*(np.cos(arg + true_anomaly) / (1 + e*np.cos(true_anomaly)))

    delta_asc = sqrt_amu * delta_v_normal*(np.sin(arg + true_anomaly) / (np.sin(i*(1 + e*np.cos(true_anomaly)))))

    arg_term1 = -delta_v_radial * np.cos(true_anomaly) / e
    arg_term2 = delta_v_angular * (2 + e * np.cos(true_anomaly) * np.sin(true_anomaly)) / (e * (1 + e * np.cos(true_anomaly)))
    arg_term3 = -delta_v_normal * np.sin(arg + true_anomaly)) / np.tan(i * (1 + e * np.cos(true_anomaly))
    delta_arg = sqrt_amu * (arg_term1 + arg_term2 + arg_term3)

    return [delta_a, delta_e, delta_i, delta_asc, delta_arg]


def get_orbital_elements(pos_vector, vel_vector, mu=398600, elements_type='standard'):
    """
    Calculates the orbital elements of an orbit based off of a body's position and velocity vectors

    Args:
        pos_vector: Position vector, a list in the form [r_x, r_y, r_z]

        vel_vector: Velocity vector, a list in the form [v_x, v_y, v_z]

        mu: Gravitational parameter, set to Earth's gravitational parameter 398,600 by default

        elements_type: Type of elements being returned, 'standard' by default, can also be set to 'curtis'
            to get orbital elements defined in Curtis' Orbital Mechanics for Engineering Students

    Returns:
        If elements_type='standard', a list of the form:
            [semi_major_axis, eccentricity, inclination, right_ascension, argument_perigee, true_anomaly]

        If elements_type='curtis', a list of the form:
            [specific_angular_momentum, inclination, right_ascension, eccentricity, argument_perigee, true_anomaly]

        All units in meters, radians, and radians/second
    """

    def magnitude(vector):
        return np.sqrt(np.sum(np.power(vector, 2)))

    distance = magnitude(pos_vector)
    speed = magnitude(vel_vector)

    h_vec = np.cross(pos_vector, vel_vector)
    specific_angular_momentum = magnitude(h_vec)

    inclination = np.arccos(h_vec[2] / specific_angular_momentum)

    N_vec = np.cross(np.array([0, 0, 1]), h_vec)
    N_mag = magnitude(N_vec)

    if N_vec[1] >= 0:
        right_ascension = np.arccos(N_vec[0] / N_mag)
    else:
        right_ascension = 2 * np.pi - np.arccos(N_vec[0] / N_mag)

    radial_velocity = np.dot(pos_vector, vel_vector) / distance
    e_vec = (1 / mu) * ((speed ** 2 - mu / distance) * pos_vector - distance * radial_velocity * vel_vector)
    eccentricity = magnitude(e_vec)

    if e_vec[2] >= 0:
        argument_perigee = np.arccos(np.dot(N_vec / N_mag, e_vec / eccentricity))
    else:
        argument_perigee = 2 * np.pi - np.arccos(np.dot(N_vec / N_mag, e_vec / eccentricity))

    if radial_velocity >= 0:
        true_anomaly = np.arccos(np.dot(e_vec / eccentricity, pos_vector / distance))
    else:
        true_anomaly = 2 * np.pi - np.arccos(np.dot(e_vec / eccentricity, pos_vector / distance))

    semi_major_axis = (distance * (1 + eccentricity * np.cos(true_anomaly))) / (1 - eccentricity ** 2)

    if elements_type == 'standard':
        return [semi_major_axis, eccentricity, inclination, right_ascension, argument_perigee, true_anomaly]
    elif elements_type == 'curtis':
        return [specific_angular_momentum, inclination, right_ascension, eccentricity, argument_perigee, true_anomaly]