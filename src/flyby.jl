"""
    flyby.jl

    This module implements gravity assist (flyby) trajectory analysis and porkchop plot computation for interplanetary missions.
    It provides data structures and functions to efficiently compute optimal transfer trajectories involving a planetary flyby, 
    including delta-v and timing analysis.

    Main features:
    - FlybyResult and FlybyWorkspace data structures for storing results and pre-allocated workspace
    - porkchop_flyby function for computing porkchop plots with a flyby maneuver
"""


using StaticArrays
using LinearAlgebra
using Base.Threads

struct FlybyResult
    tdep::Vector{Float64}
    tarr::Vector{Float64}
    best_tfly::Matrix{Float64}
    dv_total::Matrix{Float64}
    dv_dep::Matrix{Float64}
    dv_fly::Matrix{Float64}
    dv_arr::Matrix{Float64}
    c3::Matrix{Float64}
    vinf_arr::Matrix{Float64}
    ok::BitMatrix
end

# Pre-allocated workspace to prevent GC activity during hot loops
struct FlybyWorkspace
    r_dep::Vector{Vec3}
    v_dep::Vector{Vec3}
    r_fly::Vector{Vec3}
    v_fly::Vector{Vec3}
    r_arr::Vector{Vec3}
    v_arr::Vector{Vec3}
end

function FlybyWorkspace(n_dep::Int, n_fly::Int, n_arr::Int)
    FlybyWorkspace(
        Vector{Vec3}(undef, n_dep), Vector{Vec3}(undef, n_dep),
        Vector{Vec3}(undef, n_fly), Vector{Vec3}(undef, n_fly),
        Vector{Vec3}(undef, n_arr), Vector{Vec3}(undef, n_arr)
    )
end

function porkchop_flyby(sys::TwoBodySystem, 
                        body_dep::BodyModel, 
                        body_fly::BodyModel, 
                        body_arr::BodyModel,
                        tdep_range::NTuple{3,Float64},
                        tfly_range::NTuple{3,Float64},
                        tarr_range::NTuple{3,Float64};
                        R_fly_min::Float64 = body_fly.body.radius + 1000.0,
                        μ_fly::Float64 = body_fly.body.μ,
                        rpark_dep::Float64 = 6678.0,
                        rpark_arr::Float64 = NaN, 
                        tof1_min::Float64 = 50*86400.0,
                        tof1_max::Float64 = 800*86400.0,
                        tof2_min::Float64 = 50*86400.0,
                        tof2_max::Float64 = 1500*86400.0,
                        dv_cap::Float64 = 20.0)

    # Grid setup
    t0_d, t1_d, dt_d = tdep_range
    t0_f, t1_f, dt_f = tfly_range
    t0_a, t1_a, dt_a = tarr_range

    nd = Int(floor((t1_d - t0_d)/dt_d)) + 1
    nf = Int(floor((t1_f - t0_f)/dt_f)) + 1
    na = Int(floor((t1_a - t0_a)/dt_a)) + 1

    tdeps = [t0_d + (i-1)*dt_d for i in 1:nd]
    tflys = [t0_f + (i-1)*dt_f for i in 1:nf]
    tarrs = [t0_a + (i-1)*dt_a for i in 1:na]

    # Pre-allocation of output matrices
    res_dv_tot = fill(NaN, nd, na)
    res_dv_dep = fill(NaN, nd, na)
    res_dv_fly = fill(NaN, nd, na)
    res_dv_arr = fill(NaN, nd, na)
    res_c3     = fill(NaN, nd, na)
    res_vinf   = fill(NaN, nd, na)
    res_tfly   = fill(NaN, nd, na)
    res_ok     = falses(nd, na)

    # Ephemeris Pre-calculation (Single-threaded to avoid SPICE race conditions)
    ws = FlybyWorkspace(nd, nf, na)
    precompute_states!(ws.r_dep, ws.v_dep, sys, body_dep, tdeps)
    precompute_states!(ws.r_fly, ws.v_fly, sys, body_fly, tflys)
    precompute_states!(ws.r_arr, ws.v_arr, sys, body_arr, tarrs)

    μ_sun = sys.central.μ
    μ_dep = body_dep.body.μ
    μ_arr = body_arr.body.μ
    
    
    @threads :static for i in 1:nd
        r_d = ws.r_dep[i]
        v_d = ws.v_dep[i]
        t_d = tdeps[i]

        @inbounds for k in 1:na
            r_a = ws.r_arr[k]
            v_a = ws.v_arr[k]
            t_a = tarrs[k]

            best_dv_local = Inf
            
            # Inner Scan: Flyby time
            for j in 1:nf
                t_f = tflys[j]
                tof1 = t_f - t_d
                tof2 = t_a - t_f

                # Early exit on TOF constraints
                if tof1 < tof1_min || tof1 > tof1_max || tof2 < tof2_min || tof2 > tof2_max
                    continue
                end

                r_f = ws.r_fly[j]
                v_f = ws.v_fly[j]

                # Leg 1: Dep -> Flyby
                (v1_dep, v1_arr, ok1) = lambert_uv(μ_sun, r_d, r_f, tof1, false)
                if !ok1 continue end

                # Leg 2: Flyby -> Arr
                (v2_dep, v2_arr, ok2) = lambert_uv(μ_sun, r_f, r_a, tof2, false)
                if !ok2 continue end

                # Delta V Calculations
                # 1. Departure
                vinf_dep_vec = v1_dep - v_d
                vinf_dep_sq  = vinf_dep_vec ⋅ vinf_dep_vec
                vinf_dep_mag = sqrt(vinf_dep_sq)
                
                dv_launch = if isfinite(rpark_dep)
                    dv_from_vinf(μ_dep, rpark_dep, vinf_dep_mag)
                else
                    vinf_dep_mag
                end

                # 2. Flyby (Powered/Unpowered)
                vinf_in  = v1_arr - v_f
                vinf_out = v2_dep - v_f
                dv_flyby = powered_flyby_delta_v(vinf_in, vinf_out, R_fly_min, μ_fly)

                # 3. Arrival
                vinf_arr_vec = v2_arr - v_a
                vinf_arr_mag = sqrt(vinf_arr_vec ⋅ vinf_arr_vec)
                
                dv_arrival = if isfinite(rpark_arr)
                    dv_from_vinf(μ_arr, rpark_arr, vinf_arr_mag)
                else
                    vinf_arr_mag 
                end

                total_dv = dv_launch + dv_flyby + dv_arrival

                # Minimize
                if total_dv < best_dv_local
                    best_dv_local = total_dv
                    if total_dv < dv_cap
                        res_dv_tot[i,k] = total_dv
                        res_dv_dep[i,k] = dv_launch
                        res_dv_fly[i,k] = dv_flyby
                        res_dv_arr[i,k] = dv_arrival
                        res_c3[i,k]     = vinf_dep_sq
                        res_vinf[i,k]   = vinf_arr_mag
                        res_tfly[i,k]   = t_f
                        res_ok[i,k]     = true
                    end
                end
            end 
        end 
    end

    return FlybyResult(tdeps, tarrs, res_tfly, res_dv_tot, res_dv_dep, res_dv_fly, res_dv_arr, res_c3, res_vinf, res_ok)
end