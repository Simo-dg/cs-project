using StaticArrays
using LinearAlgebra

@inline function flyby_turning_angle(rp::Float64, vinf::Float64, μp::Float64)
    x = 1.0/(1.0 + rp*(vinf^2)/μp)
    x = clamp(x, 0.0, 1.0)
    return 2*asin(x)
end

function powered_flyby_delta_v(vinf_in::SVector{3,Float64}, vinf_out::SVector{3,Float64},
                               rp::Float64, μp::Float64)
    vin_in = sqrt(vinf_in⋅vinf_in)
    vin_out = sqrt(vinf_out⋅vinf_out)

    vp_in = sqrt(vin_in^2 + 2*μp/rp)
    vp_out = sqrt(vin_out^2 + 2*μp/rp)

    c = (vinf_in⋅vinf_out)/(vin_in*vin_out)
    c = clamp(c, -1.0, 1.0)
    φ = acos(c)

    return sqrt(vp_in^2 + vp_out^2 - 2*vp_in*vp_out*cos(φ))
end