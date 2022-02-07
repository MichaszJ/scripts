function hohmann_mission_profile(orbits_dict::OrderedDict, payloads_dict::Dict, vehicle_params::Dict; print_results=true, return_results=false)
    # defining constants
    μ = 398600
    r_e = 6378
    g0 = 9.80665
    
    # defining internal functions
    function _get_orbit_params(orbit_params)
        if orbit_params[1] == "C"
            return orbit_params[1], orbit_params[2], orbit_params[2]
        elseif orbit_params[1] == "E"
            return orbit_params[1], orbit_params[2], orbit_params[3]
        end
    end
    
    # defining performance parameters
    Δv_total = []
    burns = []
    payload_deployed = []
    m_fuel_used = []
    
    # vehicle parameters
    dry_mass = vehicle_params["dmass"]
    wet_mass = vehicle_params["wmass"]
    fuel_mass = wet_mass - dry_mass
    Isp = vehicle_params["Isp"]    
    total_mass = dry_mass + wet_mass + sum([i[2] for i in payload_masses])
    
    # initializing loop vars
    r_a_prev, r_p_prev, orbit_prev, type_prev = 0, 0, 0, 0
    
    # plotting function vars
    radii = []
    
    for (index, (orbit_phase, orbit_params)) in enumerate(orbits_dict)
        orbit_type, r_a, r_p = _get_orbit_params(orbit_params)
        r_a += r_e
        r_p += r_e
        Δm = 0
        Δmf = 0
        
        if index == 1
            r_a_prev, r_p_prev, orbit_prev, type_prev = r_a, r_p, orbit_phase, orbit_type
        else            
            if type_prev == "C"
                v_start = sqrt(μ / r_p_prev)
                
                if orbit_type == "C"
                    if r_a > r_a_prev
                        a_transfer = 0.5 * (r_p_prev + r_a)
                        
                        v_transfer = sqrt(2*μ) * sqrt(1/r_p_prev - 1/(2 * a_transfer))
                        v_transfer_apo = sqrt(2*μ) * sqrt(1/r_a - 1/(2 * a_transfer))
                        
                        v_circ = sqrt(μ / r_a)
                        
                        Δv_transfer = v_transfer - v_start
                        Δmf_transfer = total_mass - total_mass / exp((Δv_transfer * 1000) / (Isp * g0))  
                        push!(burns, "$(orbit_phase) transfer burn")
                        push!(Δv_total, Δv_transfer)
                        push!(m_fuel_used, Δmf_transfer)
                        push!(radii, [r_p_prev, r_a])
                        
                        total_mass -= Δmf_transfer
                        
                        Δv_circ = v_circ - v_transfer_apo
                        Δmf_circ = total_mass - total_mass / exp((Δv_circ * 1000) / (Isp * g0))  
                        push!(burns, "$(orbit_phase) circularization burn")
                        push!(Δv_total, Δv_circ)
                        push!(m_fuel_used, Δmf_circ)
                        push!(radii, [r_p, r_a])
                        
                        total_mass -= Δmf_circ
                    else
                        a_transfer = 0.5 * (r_p + r_a_prev)
                        
                        v_transfer = sqrt(2*μ) * sqrt(1/r_a_prev - 1/(2 * a_transfer))
                        v_transfer_per = sqrt(2*μ) * sqrt(1/r_p - 1/(2 * a_transfer))
                        
                        v_circ = sqrt(μ / r_p)
                        
                        Δv_transfer = v_start - v_transfer
                        Δmf_transfer = total_mass - total_mass / exp((Δv_transfer * 1000) / (Isp * g0))  
                        push!(burns, "$(orbit_phase) transfer burn")
                        push!(Δv_total, Δv_transfer)
                        push!(m_fuel_used, Δmf_transfer)
                        push!(radii, [r_a_prev, r_p])
                        
                        total_mass -= Δmf_transfer
                        
                        Δv_circ = v_transfer_per - v_circ
                        Δmf_circ = total_mass - total_mass / exp((Δv_circ * 1000) / (Isp * g0))  
                        push!(burns, "$(orbit_phase) circularization burn")
                        push!(Δv_total, Δv_circ)
                        push!(Δv_total, Δmf_circ)
                        push!(radii, [r_p, r_a])
                        
                        total_mass -= Δmf_circ
                    end

                elseif orbit_type == "E"                    
                    # establish apoapsis
                    if r_a != r_a_prev
                        a_transfer = 0.5 * (r_p_prev + r_a)
                        v_transfer = sqrt(2*μ) * sqrt(1/r_p_prev - 1/(2 * a_transfer))
                        
                        Δv_transfer = abs(v_transfer - v_start)
                        Δmf_transfer = total_mass - total_mass / exp((Δv_transfer * 1000) / (Isp * g0))  
                        push!(burns, "$(orbit_phase) apoapsis establishment burn")
                        push!(Δv_total, Δv_transfer)
                        push!(m_fuel_used, Δmf_transfer)
                        push!(radii, [r_p_prev, r_a])
                        
                        total_mass -= Δmf_transfer
                    end
                    
                    # establish periapsis
                    if r_p != r_p_prev
                        a_transfer = 0.5 * (r_p_prev + r_a)
                        v_apo = sqrt(2*μ) * sqrt(1/r_a - 1/(2 * a_transfer))
                        
                        a_orbit = 0.5 * (r_p + r_a)
                        v_transfer = sqrt(2*μ) * sqrt(1/r_a - 1/(2 * a_orbit))
                        
                        Δv_transfer = abs(v_transfer - v_apo)
                        Δmf_transfer = total_mass - total_mass / exp((Δv_transfer * 1000) / (Isp * g0))  
                        push!(burns, "$(orbit_phase) periapsis establishment burn")
                        push!(Δv_total, Δv_transfer)
                        push!(m_fuel_used, Δmf_transfer)
                        push!(radii, [r_p, r_a])
                        
                        total_mass -= Δmf_transfer
                    end
                end
            
            elseif type_prev == "E"
                if orbit_type == "C"
                    # establish apoapsis
                    if r_a != r_a_prev
                        a_prev = 0.5 * (r_p_prev + r_a_prev)
                        v_per = sqrt(2*μ) * sqrt(1/r_p_prev - 1/(2 * a_prev))

                        a_transfer = 0.5 * (r_p_prev + r_a)
                        v_transfer = sqrt(2*μ) * sqrt(1/r_p - 1/(2 * a_transfer))

                        Δv_transfer = abs(v_transfer - v_per)
                        Δmf_transfer = total_mass - total_mass / exp((Δv_transfer * 1000) / (Isp * g0))
                        push!(burns, "$(orbit_phase) transfer burn")
                        push!(Δv_total, Δv_transfer)
                        push!(m_fuel_used, Δmf_transfer)
                        push!(radii, [r_p, r_a])

                        total_mass -= Δmf_transfer
                    end

                    # circularize
                    v_apo = sqrt(2*μ) * sqrt(1/r_a - 1/(2 * a_transfer))
                    v_circ = sqrt(μ / r_a)

                    Δv_circ = abs(v_circ - v_apo)
                    Δmf_circ = total_mass - total_mass / exp((Δv_circ * 1000) / (Isp * g0))
                    push!(burns, "$(orbit_phase) circularization burn")
                    push!(Δv_total, Δv_circ)
                    push!(m_fuel_used, Δmf_transfer)
                    push!(radii, [r_p, r_a])

                    total_mass -= Δmf_transfer
                    
                elseif orbit_type == "E"
                    # establish apoapsis
                    if r_a != r_a_prev
                        a_prev = 0.5 * (r_p_prev + r_a_prev)
                        v_per = sqrt(2*μ) * sqrt(1/r_p_prev - 1/(2 * a_prev))

                        a_transfer = 0.5 * (r_p_prev + r_a)
                        v_transfer = sqrt(2*μ) * sqrt(1/r_p - 1/(2 * a_transfer))

                        Δv_transfer = abs(v_transfer - v_per)
                        Δmf_transfer = total_mass - total_mass / exp((Δv_transfer * 1000) / (Isp * g0))
                        push!(burns, "$(orbit_phase) transfer burn")
                        push!(Δv_total, Δv_transfer)
                        push!(m_fuel_used, Δmf_transfer)
                        push!(radii, [r_p, r_a])

                        total_mass -= Δmf_transfer
                    end
                    
                    # establish periapsis
                    if r_p != r_p_prev
                        a_transfer = 0.5 * (r_p_prev + r_a)
                        v_apo = sqrt(2*μ) * sqrt(1/r_a - 1/(2 * a_transfer))
                        
                        a_new = 0.5 * (r_a + r_p)
                        v_per = sqrt(2*μ) * sqrt(1/r_a - 1/(2 * a_new))
                        
                        Δv_transfer = abs(v_apo - v_per)
                        Δmf_transfer = total_mass - total_mass / exp((Δv_transfer * 1000) / (Isp * g0))
                        push!(burns, "$(orbit_phase) periapsis establishment burn")
                        push!(Δv_total, Δv_transfer)
                        push!(m_fuel_used, Δmf_transfer)
                        push!(radii, [r_p, r_a])

                        total_mass -= Δmf_transfer
                    end
                end
            end
        end
        
        if haskey(payloads_dict, orbit_phase)
            total_mass -= payloads_dict[orbit_phase]
        end
        
        r_a_prev, r_p_prev, orbit_prev, type_prev = r_a, r_p, orbit_phase, orbit_type
    end
    
    if print_results
        for (burn, Δv, Δmf) in zip(burns, Δv_total, m_fuel_used)
            println("$burn")
            println("\tΔv = $(round(Δv * 1000, digits=2)) m/s")
            println("\tFuel used = $(round(Δmf, digits=2)) kg\n")
        end

        total_fuel = sum(m_fuel_used)

        println("Total Δv = $(round(sum(Δv_total), digits=3) * 1000) m/s")
        println("Fuel required = $((round(total_fuel, digits=3))) kg")
        println("Fuel remaining = $(round(fuel_mass - total_fuel, digits=2)) / $fuel_mass kg")
    end
    
    if return_results
        return burns, Δv_total, radii
    end
end