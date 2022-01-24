using LinearAlgebra, Dates
include("orbital_dynamics.jl")

function determine_eccentric_anomaly(eccentricity, mean_anomaly; tolerance=1e-8)
    if mean_anomaly < pi
        eccentric_anomaly = mean_anomaly + eccentricity/2
    else
        eccentric_anomaly = mean_anomaly - eccentricity/2
    end

    function ratio(eccentricty, mean_anomaly, eccentric_anomaly)
        return (eccentric_anomaly - eccentricity*sin(eccentric_anomaly) - mean_anomaly) / (1 - eccentricity*cos(eccentric_anomaly))
    end

    current_ratio = ratio(eccentricity, mean_anomaly, eccentric_anomaly)

    while abs(current_ratio) > tolerance
        eccentric_anomaly = eccentric_anomaly - current_ratio
        current_ratio = ratio(eccentricity, mean_anomaly, eccentric_anomaly)
    end

    return eccentric_anomaly
end

function geocentric_state_vector_transform(angular_momentum, orbital_elements; mu=398600)
    h = angular_momentum
    e, i, asc, arg, theta = orbital_elements

    p⃗_perifocal = ((h^2 / mu) / (1 + e * cos(theta))) * [cos(theta); sin(theta); 0]
    v⃗_perifocal = (mu / h) * [-sin(theta); e + cos(theta); 0]

    transform = [
        [-sin(asc) * cos(i) * sin(arg) + cos(asc) * cos(arg)] [-sin(asc) * cos(i) * cos(arg) - cos(asc) * sin(arg)] [sin(asc) * sin(i)]
        [cos(asc) * cos(i) * sin(arg) + sin(asc) * cos(arg)] [cos(asc) * cos(i) * cos(arg) - sin(asc) * sin(arg)] [-cos(asc) * sin(i)]
        [sin(i) * sin(arg)] [sin(i) * cos(arg)] [cos(i)]
    ]

    pos_vec_geocentric = transform * p⃗_perifocal
    vel_vec_geocentric = transform * v⃗_perifocal

    return pos_vec_geocentric, vel_vec_geocentric
end

function position_to_asc_dec(position_vector)
    magnitude = sqrt(sum(position_vector.^2))

    direct_l = position_vector[1] / magnitude
    direct_m = position_vector[2] / magnitude
    direct_n = position_vector[3] / magnitude

    declination = asin(direct_n)

    direct_m > 0 ? ascension = acos(direct_l / cos(declination)) : ascension = 2*pi - acos(direct_l / cos(declination))

    return ascension, declination
end

function ground_track(orbital_elements, r_apo, r_per, t_init, t_final; num_steps=100, return_times=false, mu=398600, radius=6378, j2=1.08263e-3, body_angular_vel=7.292124e-5)
    a, e, i, asc, arg, theta = orbital_elements

    common_term = (3 / 2) * ((sqrt(mu) * j2 * radius^2) / ((1 - e^2) * a^(7/2)))
    delta_right_ascension = -common_term * cos(i)
    delta_argument_perigee = -common_term * ((5 / 2) * sin(i)^2 - 2)

    angular_momentum = sqrt(mu * r_per * (1 + e))
    period = (a^(3/2)) * 2*pi / sqrt(mu)

    eccentric_anomaly_0 = 2 * atan(tan(theta / 2) * sqrt((1 - e) / (1 + e)))
    mean_anomaly_0 = eccentric_anomaly_0 - e * sin(eccentric_anomaly_0)
    t_0 = period * mean_anomaly_0 / (2 * pi)

    longitude_vec = []
    latitude_vec = []
    times = []

    earth_rot = 2 * pi / (24 * 60 * 60)

    step_size = (t_final - t_init) / num_steps
    for c=1:num_steps
        t_i = t_0 + c * step_size

        mean_anomaly_i = t_i * (2 * pi) / period
        eccentric_anomaly_i = determine_eccentric_anomaly(e, mean_anomaly_i)
        true_anomaly_i = 2 * atan(tan(eccentric_anomaly_i / 2) * sqrt((1 + e) / (1 - e)))

        right_ascension_i = asc + step_size * delta_right_ascension
        argument_perigee_i = arg + step_size * delta_argument_perigee

        pos_geo, vel_geo = geocentric_state_vector_transform(angular_momentum, [e, i, right_ascension_i, argument_perigee_i, true_anomaly_i])

        transform_theta = step_size * body_angular_vel
        rotation_matrix = [
            [cos(transform_theta)] [sin(transform_theta)] [0]
            [-sin(transform_theta)] [cos(transform_theta)] [0]
            [0] [0] [1]
        ]

        pos_geo_rotating = rotation_matrix * pos_geo[:,:]

        longitude, latitude = position_to_asc_dec(pos_geo_rotating)

        push!(longitude_vec, longitude - (c * step_size) * earth_rot)
        push!(latitude_vec, latitude)
        push!(times, t_i)

    end

    if return_times 
        return longitude_vec, latitude_vec, times 
    else 
        return longitude_vec, latitude_vec
    end
end

function lamberts_problem(r⃗₁, r⃗₂, Δt; trajectory="prograde", μ=398600, z_init=1, rel_error=1e-3)
    @assert trajectory == "prograde" || trajectory == "retrograde" "Trajectory should be either \"prograde\" or \"trajectory\""
    
    r₁ = sqrt(r⃗₁ ⋅ r⃗₁)
    r₂ = sqrt(r⃗₂ ⋅ r⃗₂)

    r_cross = r⃗₁ × r⃗₂

    if trajectory == "prograde"
        if r_cross[3] ≥ 0
            Δθ = acos((r⃗₁ ⋅ r⃗₂) / (r₁ * r₂))
        else
            Δθ = 2*π - acos((r⃗₁ ⋅ r⃗₂) / (r₁ * r₂))
        end
    elseif trajectory == "retrograde"
        if r_cross[3] ≥ 0 
            Δθ = 2*π - acos((r⃗₁ ⋅ r⃗₂) / (r₁ * r₂))
        else
            Δθ = acos((r⃗₁ ⋅ r⃗₂) / (r₁ * r₂))
        end
    end

    A = sin(Δθ) * sqrt(r₁ * r₂ / (1 - cos(Δθ)))

    function _stumpff_S(z)
        if z > 0
            S = (√z - sin(√z)) / (√z ^ 3)
        elseif z < 0
            S = (sinh(√-z) - √-z) / (√-z ^ 3)
        else
            S = 1/6
        end
    
        return S
    end
    
    function _stumpff_C(z)
        if z > 0
            C = (1 - cos(√z)) / z
        elseif z < 0
            C = (cosh(√-z) - 1) / -z
        else
            C = 1/2
        end
    
        return C
    end
    
    _stumpff_S_prime(z) = (1 / (2*z)) * (_stumpff_C(z) - 3 * _stumpff_S(z))
    _stumpff_C_prime(z) = (1 / (2*z)) * (1 - z * _stumpff_S(z) - 2 * _stumpff_C(z))

    y(z) = r₁ + r₂ + A * (z * _stumpff_S(z) - 1) /  √_stumpff_C(zᵢ)
    y_prime(z) = (A / 4) * √_stumpff_C(z)
    
    F(z) = ((y(z) / _stumpff_C(z))^(3/2)) * _stumpff_S(z) + A * √y(z) - (√μ) * Δt
    
    function F_prime(z)
        if z != 0
            term1 = (y(z) / _stumpff_C(z))^(3/2)
            term2 = ((_stumpff_C(z) - (3/2) * (_stumpff_S(z) / _stumpff_C(z))) + (3/4) * (_stumpff_S(z)^2 / _stumpff_C(z))) / (2*z)
            term3 = (A / 8) * (3 * (_stumpff_S(z) / _stumpff_C(z)) * √y(z) + A * √(_stumpff_C(z) / y(z)))
            F′ = term1 * term2 + term3
        else
            F′ = (√2 / 40) * y(0)^(3/2) + (A / 8) * (√y(0) + A * √(1 / (2 * y(0))))
        end

        return F′
    end

    err = 100
    zᵢ = z_init
    while err > rel_error
        z_prev = zᵢ
        zᵢ = zᵢ - F(zᵢ) / F_prime(zᵢ)

        err = abs((zᵢ - z_prev) / z_prev)
    end

    yᵢ = y(zᵢ)

    lagrange_f = 1 - yᵢ / r₁
    lagrange_g = A * √(yᵢ / μ)
    lagrange_ḟ = (√μ / (r₁ * r₂)) * √(yᵢ / _stumpff_C(zᵢ)) * (zᵢ * _stumpff_S(zᵢ) - 1)
    lagrange_ġ = 1 - yᵢ / r₂

    v⃗₁ = (r⃗₂ - lagrange_f .* r⃗₁) ./ lagrange_g
    v⃗₂ = (lagrange_ġ .* r⃗₂ - r⃗₁) ./ lagrange_g

    return v⃗₁, v⃗₂
end

function get_planet_ephemeris(datetime; planet="Earth", μ=1.327e11)
    #year, month, day = date
    #hour, minute, second = time
    year, month, day = Dates.yearmonthday(datetime)
    hour = Dates.hour(datetime)
    minute = Dates.minute(datetime)
    second = Dates.second(datetime)

    fterm1 = convert(Int64, floor((month + 9)/12))
    fterm2 = convert(Int64, floor((7 * (year + fterm1)) / 4))
    fterm3 = convert(Int64, floor(275*month / 9))
    J₀ = 367*year - fterm2 + fterm3 + day + 1721013.5

    UT = hour + minute/60 + second/3600

    JD = J₀ + UT/24

    T₀ = (JD - 2451545) / 36525

    planet_elements = Dict(
        "Mercury" => [[0.38709927, 0.00000037], [0.20563593, 0.00001906], [7.00497902, -0.00594749], [48.33076593, -0.12534081], [77.45779628, 0.16047689], [252.25032350, 149472.67411175]],
        "Venus" => [[0.72333566, 0.00000390], [0.00677672, -0.00004107], [3.39467605, -0.00078890], [76.67984255, -0.27769418], [131.60246718, 0.00268329], [181.97909950, 58517.81538729]],
        "Earth" => [[1.00000261, 0.00000562], [0.01671123, -0.00004392], [-0.00001531, -0.01294668], [0.0, 0.0], [102.93768193, 0.32327364], [100.46457166, 35999.37244981]],
        "Mars" => [[1.52371034, 0.0001847], [0.09339410, 0.00007882], [1.84969142, -0.00813131], [49.55953891, -0.29257343], [-23.94362959, 0.44441088], [-4.55343205, 19140.30268499]],
        "Jupiter" => [[5.20288700, -0.00011607], [0.04838624, -0.00013253], [1.30439695, -0.00183714], [100.47390909, 0.20469106], [14.72847983, 0.21252668], [34.39644501, 3034.74612775]],
        "Saturn" => [[9.53667594, -0.00125060], [0.05386179, -0.00050991], [2.48599187, 0.00193609], [113.66242448, -0.28867794], [92.59887831, -0.41897216], [49.95424423, 1222.49362201]],
        "Uranus" => [[19.18916464, -0.00196176], [0.04725744, -0.00004397], [0.77263783, -0.00242939], [74.01692503, 0.04240589], [170.95427630, 0.40805281], [313.23810451, 428.48202785]],
        "Neptune" => [[30.06992276, 0.00026291], [0.00859048, 0.00005105], [1.77004347, 0.00035372], [131.78422574, -0.00508664], [44.96476227, -0.32241464], [-55.12002969, 218.45945325]],
        "Pluto" => [[39.48211675, -0.00031596], [0.24882730, 0.00005170], [17.14001206, 0.00004818], [110.30393684, -0.01183482], [224.06891629, -0.04062942], [238.92903833, 145.20780515]]
    )

    elements_init = planet_elements[planet]
    elements = [element[1] + element[2] * T₀ for element in elements_init]
    
    # converting degree units to radians and wrapping angles
    function _wrap_angle(angle; type="rad")
        @assert type == "rad" || type == "deg" "Invalid type $(type), choose either rad or deg"
        
        new_angle = copy(angle)
        if type == "rad"
            if new_angle > 0
                while new_angle > 2*π
                    new_angle -= 2*π
                end
            else
                while new_angle < 0
                    new_angle += 2*π
                end
            end
        else
            if new_angle > 0
                while new_angle > 360
                    new_angle -= 360
                end
            else
                while new_angle < 0
                    new_angle += 360
                end
            end
        end

        return new_angle
    end

    elements[3:6] = [_wrap_angle(deg2rad(elements[i])) for i ∈ 3:6]

    AU = 149597870.7
    elements[1] = elements[1] * AU

    a, e, i, Ω, ω̄, L = elements

    h = sqrt(μ * a * (1 - e^2))
    ω = _wrap_angle(ω̄ - Ω)
    M = _wrap_angle(L - ω̄)

    E = _wrap_angle(determine_eccentric_anomaly(e, M))

    θ = _wrap_angle(2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)))

    pos_vec_geocentric, vel_vec_geocentric = geocentric_state_vector_transform(h, [e, i, Ω, ω, θ], mu = μ)

    return pos_vec_geocentric, vel_vec_geocentric
end

function interplanetary_trajectory(departure_datetime, arrival_datetime, departure_planet, arrival_planet)
    pos_departure, vel_departure = get_planet_ephemeris(departure_datetime, planet=departure_planet)
    pos_arrival, vel_arrival = get_planet_ephemeris(arrival_datetime, planet=arrival_planet)

    time_duration = convert(Dates.Second, arrival_datetime - departure_datetime).value
    departure_velocity, arrival_velocity = lamberts_problem(pos_departure, pos_arrival, time_duration, μ=1.327e11)

    transfer_elements = get_orbital_elements(pos_departure, departure_velocity, μ=1.327e11)

    return transfer_elements
end