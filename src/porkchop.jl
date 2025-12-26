# src/porkchop.jl
using StaticArrays
using LinearAlgebra

struct PorkchopResult
    tdep::Vector{Float64}          # seconds
    tarr::Vector{Float64}          # seconds
    dv::Matrix{Float64}            # total Δv (km/s) OR proxy (km/s)
    c3::Matrix{Float64}            # departure C3 (km^2/s^2)
    vinf_arr::Matrix{Float64}      # arrival v∞ (km/s)
    ok::BitMatrix
end

# Δv from circular parking orbit radius rpark about planet μp with given v∞
@inline function dv_from_vinf(μp::Float64, rpark::Float64, vinf::Float64)
    vcirc = sqrt(μp / rpark)
    vesc  = sqrt(2*μp / rpark)
    return sqrt(vinf*vinf + vesc*vesc) - vcirc
end

# Precompute ephemeris states for a vector of times.
# CRITICAL: if the ephemeris is SPICE-based, this keeps all SPICE calls serial
# and the threaded porkchop scan becomes safe + much faster.
function precompute_states(sys::TwoBodySystem, bm::BodyModel, ts::Vector{Float64})
    rs = Vector{Vec3}(undef, length(ts))
    vs = Vector{Vec3}(undef, length(ts))
    @inbounds for k in eachindex(ts)
        rs[k], vs[k] = state(sys, bm, ts[k])
    end
    return rs, vs
end

"""
    porkchop_grid_dep_arr(sys, dep, arr, tdep_range, tarr_range;
                          longway=false,
                          metric=:dv,
                          dep_μ=dep.body.μ, dep_rpark=NaN,
                          arr_μ=arr.body.μ, arr_rpark=NaN,
                          tof_min=60*86400.0, tof_max=500*86400.0,
                          dv_cap=Inf)

TRUE porkchop grid:
- x-axis: departure time samples
- y-axis: arrival time samples
- value:
    metric=:dv    -> Δv_departure(from parking) + Δv_arrival(to capture)
                     (requires dep_rpark and arr_rpark finite; otherwise falls back to proxy)
    metric=:vinf  -> v∞dep + v∞arr proxy (km/s)

Threaded over departure index. SPICE-safe because ephemerides are precomputed serially.
"""
function porkchop_grid_dep_arr(sys::TwoBodySystem, dep::BodyModel, arr::BodyModel,
                               tdep_range::NTuple{3,Float64},
                               tarr_range::NTuple{3,Float64};
                               longway::Bool=false,
                               metric::Symbol=:dv,
                               dep_μ::Float64=dep.body.μ,
                               dep_rpark::Float64=NaN,
                               arr_μ::Float64=arr.body.μ,
                               arr_rpark::Float64=NaN,
                               tof_min::Float64=60*86400.0,
                               tof_max::Float64=500*86400.0,
                               dv_cap::Float64=Inf)

    t0, t1, dt = tdep_range
    a0, a1, da = tarr_range

    nt = Int(floor((t1 - t0)/dt)) + 1
    na = Int(floor((a1 - a0)/da)) + 1

    tdep = Vector{Float64}(undef, nt)
    tarr = Vector{Float64}(undef, na)
    @inbounds for i in 1:nt
        tdep[i] = t0 + (i-1)*dt
    end
    @inbounds for j in 1:na
        tarr[j] = a0 + (j-1)*da
    end

    dv = Matrix{Float64}(undef, nt, na)
    c3 = Matrix{Float64}(undef, nt, na)
    vinf_arr = Matrix{Float64}(undef, nt, na)
    ok = BitMatrix(undef, nt, na)

    μc = sys.central.μ

    use_parking = (metric === :dv) && isfinite(dep_rpark) && isfinite(arr_rpark)

    # --- PRECOMPUTE ephemeris states (serial -> SPICE-safe, and faster)
    r1s, v1ps = precompute_states(sys, dep, tdep)
    r2s, v2ps = precompute_states(sys, arr, tarr)

    Threads.@threads for i in 1:nt
        td = tdep[i]
        r1  = r1s[i]
        v1p = v1ps[i]

        @inbounds for j in 1:na
            ta = tarr[j]
            tof = ta - td

            # enforce TOF window (this is what makes a readable porkchop and avoids insane dv)
            if tof < tof_min || tof > tof_max
                ok[i,j] = false
                dv[i,j] = NaN
                c3[i,j] = NaN
                vinf_arr[i,j] = NaN
                continue
            end

            r2  = r2s[j]
            v2p = v2ps[j]

            v1t, v2t, success = lambert_uv(μc, r1, r2, tof; longway=longway)
            ok[i,j] = success

            if !success
                dv[i,j] = NaN
                c3[i,j] = NaN
                vinf_arr[i,j] = NaN
                continue
            end

            vdep_inf_vec = v1t - v1p
            varr_inf_vec = v2t - v2p

            c3ij = vdep_inf_vec ⋅ vdep_inf_vec
            vinf_dep = sqrt(c3ij)
            vinfA = sqrt(varr_inf_vec ⋅ varr_inf_vec)

            c3[i,j] = c3ij
            vinf_arr[i,j] = vinfA

            # metric
            if metric === :vinf
                dvij = vinf_dep + vinfA
            elseif use_parking
                dvdep = dv_from_vinf(dep_μ, dep_rpark, vinf_dep)
                dvarr = dv_from_vinf(arr_μ, arr_rpark, vinfA)
                dvij = dvdep + dvarr
            else
                dvij = vinf_dep + vinfA
            end

            # cap huge dv (mostly for visualization / masking)
            if dvij > dv_cap
                ok[i,j] = false
                dv[i,j] = NaN
            else
                dv[i,j] = dvij
            end
        end
    end

    return PorkchopResult(tdep, tarr, dv, c3, vinf_arr, ok)
end