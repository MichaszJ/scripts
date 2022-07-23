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
    arg_term3 = -delta_v_normal * np.sin(arg + true_anomaly) / np.tan(i * (1 + e * np.cos(true_anomaly)))
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


def determine_eccentric_anomaly(eccentricity, mean_anomaly, tolerance=1e-8):
    if mean_anomaly < np.pi:
        eccentric_anomaly = mean_anomaly + eccentricity/2
    else:
        eccentric_anomaly = mean_anomaly - eccentricity/2

    def ratio(eccentricty, mean_anomaly, eccentric_anomaly):
        return (eccentric_anomaly - eccentricity*np.sin(eccentric_anomaly) - mean_anomaly) / (1 - eccentricity*np.cos(eccentric_anomaly))

    current_ratio = ratio(eccentricity, mean_anomaly, eccentric_anomaly)

    while abs(current_ratio) > tolerance:
        eccentric_anomaly = eccentric_anomaly - current_ratio
        current_ratio = ratio(eccentricity, mean_anomaly, eccentric_anomaly)

    return eccentric_anomaly


def geocentric_state_vector_transform(angular_momentum, orbital_elements, mu=398600):
    h = angular_momentum
    e, i, asc, arg, theta = orbital_elements

    pos_vec_perifocal = ((h ** 2 / mu) / (1 + e * np.cos(theta))) * np.array([[np.cos(theta)], [np.sin(theta)], [0]])
    vel_vec_perifocal = (mu / h) * np.array([[-np.sin(theta)], [e + np.cos(theta)], [0]])

    transform_1 = [
        -np.sin(asc) * np.cos(i) * np.sin(arg) + np.cos(asc) * np.cos(arg),
        -np.sin(asc) * np.cos(i) * np.cos(arg) - np.cos(asc) * np.sin(arg),
        np.sin(asc) * np.sin(i)
    ]

    transform_2 = [
        np.cos(asc) * np.cos(i) * np.sin(arg) + np.sin(asc) * np.cos(arg),
        np.cos(asc) * np.cos(i) * np.cos(arg) - np.sin(asc) * np.sin(arg),
        -np.cos(asc) * np.sin(i)
    ]

    transform_3 = [
        np.sin(i) * np.sin(arg),
        np.sin(i) * np.cos(arg),
        np.cos(i)
    ]

    transform = np.array([transform_1, transform_2, transform_3])

    pos_vec_geocentric = np.matmul(transform, pos_vec_perifocal)
    vel_vec_geocentric = np.matmul(transform, vel_vec_perifocal)

    return pos_vec_geocentric, vel_vec_geocentric


def position_to_asc_dec(position_vector):
    magnitude = np.sqrt(np.sum(np.power(position_vector, 2)))

    direct_l = position_vector[0] / magnitude
    direct_m = position_vector[1] / magnitude
    direct_n = position_vector[2] / magnitude

    declination = np.arcsin(direct_n)

    if direct_m > 0:
        ascension = np.arccos(direct_l / np.cos(declination))
    else:
        ascension = 2 * np.pi - np.arccos(direct_l / np.cos(declination))

    return ascension, declination


def ground_track(orbital_elements,
                 r_apo,
                 r_per,
                 t_init,
                 t_final,
                 num_steps=100,
                 return_times=False,
                 mu=398600,
                 radius=6378,
                 j2=1.08263e-3,
                 body_angular_vel=7.292124e-5):
    """
    Returns the longitude and latitude (in radians) ground track of a satellite based on input orbital elements
    """

    a, e, i, asc, arg, theta = orbital_elements

    common_term = (3 / 2) * ((np.sqrt(mu) * j2 * radius ** 2) / ((1 - e ** 2) * np.power(a, 7 / 2)))
    delta_right_ascension = -common_term * np.cos(i)
    delta_argument_perigee = -common_term * ((5 / 2) * np.sin(i) ** 2 - 2)

    angular_momentum = np.sqrt(mu * r_per * (1 + e))
    period = np.power(a, 3 / 2) * (2 * np.pi) / np.sqrt(mu)

    eccentric_anomaly_0 = 2 * np.arctan(np.tan(theta / 2) * np.sqrt((1 - e) / (1 + e)))
    mean_anomaly_0 = eccentric_anomaly_0 - e * np.sin(eccentric_anomaly_0)
    t_0 = period * mean_anomaly_0 / (2 * np.pi)

    longitude_vec = []
    latitude_vec = []
    times = []

    earth_rot = 2 * np.pi / (24 * 60 * 60)

    step_size = (t_final - t_init) / num_steps
    for c in range(num_steps):
        t_i = t_0 + c * step_size

        mean_anomaly_i = t_i * (2 * np.pi) / period
        eccentric_anomaly_i = determine_eccentric_anomaly(e, mean_anomaly_i)
        true_anomaly_i = 2 * np.arctan(np.tan(eccentric_anomaly_i / 2) * np.sqrt((1 + e) / (1 - e)))

        right_ascension_i = asc + step_size * delta_right_ascension
        argument_perigee_i = arg + step_size * delta_argument_perigee

        pos_geo, vel_geo = geocentric_state_vector_transform(angular_momentum,
                                                             [e, i, right_ascension_i, argument_perigee_i,
                                                              true_anomaly_i])

        transform_theta = step_size * body_angular_vel
        rotation_matrix = np.array([
            [np.cos(transform_theta), np.sin(transform_theta), 0],
            [-np.sin(transform_theta), np.cos(transform_theta), 0],
            [0, 0, 1]
        ])

        pos_geo_rotating = np.matmul(rotation_matrix, pos_geo)

        longitude, latitude = position_to_asc_dec(pos_geo_rotating)

        longitude_vec.append(longitude[0] - (c * step_size) * earth_rot)
        latitude_vec.append(latitude[0])
        times.append(t_i)

    if return_times:
        return np.array(longitude_vec), np.array(latitude_vec), np.array(times)

    else:
        return np.array(longitude_vec), np.array(latitude_vec)