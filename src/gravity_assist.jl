using StaticArrays
using LinearAlgebra

@inline function flyby_turning_angle(rp::Float64, vinf::Float64, μp::Float64)
    e = 1.0 + (rp * vinf^2) / μp
    return 2.0 * asin(1.0 / e)
end

@inline function powered_flyby_delta_v(vinf_in::SVector{3,Float64}, 
                                       vinf_out::SVector{3,Float64},
                                       rp_min::Float64, 
                                       μp::Float64)
    vin_sq  = vinf_in ⋅ vinf_in
    vout_sq = vinf_out ⋅ vinf_out
    vin     = sqrt(vin_sq)
    vout    = sqrt(vout_sq)

    # Required turning angle
    cos_phi = clamp((vinf_in ⋅ vinf_out) / (vin * vout), -1.0, 1.0)
    phi_req = acos(cos_phi)

    # Max passive turning angle 
    term_in  = 1.0 + (rp_min * vin_sq) / μp
    term_out = 1.0 + (rp_min * vout_sq) / μp
    delta_max = asin(1.0/term_in) + asin(1.0/term_out)

    if phi_req <= delta_max
        # Gravity is sufficient
        return abs(vout - vin)
    else
        # Gravity insufficientx)
        vp_in  = sqrt(vin_sq + 2*μp/rp_min)
        vp_out = sqrt(vout_sq + 2*μp/rp_min)
        
        # Vector cosine law at periapsis for the unachievable angle portion
        excess_angle = phi_req - delta_max
        return sqrt(vp_in^2 + vp_out^2 - 2*vp_in*vp_out*cos(excess_angle))
    end
end