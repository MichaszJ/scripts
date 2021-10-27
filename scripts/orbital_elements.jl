using LinearAlgebra

function get_orbital_elements(pos_vector, vel_vector; mu=398600, elements_type="standard")
    function magnitude(vector)
        return sqrt(sum(vector .^ 2))
    end

    distance = magnitude(pos_vector)
    speed = magnitude(vel_vector)

    h_vec = cross(pos_vector, vel_vector)
    specific_angular_momentum = magnitude(h_vec)

    inclination = acos(h_vec[3] / specific_angular_momentum)

    N_vec = cross([0, 0, 1], h_vec)
    N_mag = magnitude(N_vec)

    N_vec[2] >= 0 ? right_ascension = acos(N_vec[1] / N_mag) : right_ascension = 2 * pi - acos(N_vec[1] / N_mag)

    radial_velocity = dot(pos_vector, vel_vector) / distance
    e_vec = (1 / mu) * ((speed^2 - mu / distance) * pos_vector - distance * radial_velocity * vel_vector)
    eccentricity = magnitude(e_vec)

    e_vec[3] >= 0 ? argument_perigee = acos(dot(N_vec / N_mag, e_vec / eccentricity)) : argument_perigee = 2 * np.pi - acos(dot(N_vec / N_mag, e_vec / eccentricity))

    radial_velocity >= 0 ? true_anomaly = acos(dot(e_vec / eccentricity, pos_vector / distance)) : true_anomaly = 2 * pi - acos(dot(e_vec / eccentricity, pos_vector / distance))

    semi_major_axis = (distance * (1 + eccentricity * cos(true_anomaly))) / (1 - eccentricity^2)

    if elements_type == "standard"
        return [semi_major_axis, eccentricity, inclination, right_ascension, argument_perigee, true_anomaly]
    elseif elements_type == "curtis"
        return [specific_angular_momentum, inclination, right_ascension, eccentricity, argument_perigee, true_anomaly]
    end
end