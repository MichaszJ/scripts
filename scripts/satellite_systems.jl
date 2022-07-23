using LinearAlgebra

function panel_in_sun(r⃗_sat, r̂_sun, r̂_panel)
	in_sun = false

	if dot(r⃗_sat, r̂_sun) < 0
		# check if spacecraft is outside eclipse region
		if magnitude(r⃗_sat - (dot(r⃗_sat, r̂_sun)) .* r̂_sun) > 6378.0
			# check if panel has sun shining on it
			if vec_angle(r̂_panel, r⃗_sun) < π/2
				in_sun = true
			end
		end
	else
		# check if panel has sun shining on it
		if vec_angle(r̂_panel, r̂_sun) < π/2
			in_sun = true
		end
	end

	return in_sun
end