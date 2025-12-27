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

mutable struct PorkchopWorkspace
    tdep::Vector{Float64}
    tarr::Vector{Float64}
    dv::Matrix{Float64}
    c3::Matrix{Float64}
    vinf_arr::Matrix{Float64}
    ok::BitMatrix
    r1s::Vector{Vec3}
    v1ps::Vector{Vec3}
    r2s::Vector{Vec3}
    v2ps::Vector{Vec3}
end

function PorkchopWorkspace(nt::Int, na::Int)
    PorkchopWorkspace(
        Vector{Float64}(undef, nt),
        Vector{Float64}(undef, na),
        Matrix{Float64}(undef, nt, na),
        Matrix{Float64}(undef, nt, na),
        Matrix{Float64}(undef, nt, na),
        BitMatrix(undef, nt, na),
        Vector{Vec3}(undef, nt),
        Vector{Vec3}(undef, nt),
        Vector{Vec3}(undef, na),
        Vector{Vec3}(undef, na),
    )
end

function ensure_workspace!(ws::PorkchopWorkspace, nt::Int, na::Int)
    if length(ws.tdep) != nt
        resize!(ws.tdep, nt)
        resize!(ws.r1s, nt)
        resize!(ws.v1ps, nt)
    end
    if length(ws.tarr) != na
        resize!(ws.tarr, na)
        resize!(ws.r2s, na)
        resize!(ws.v2ps, na)
    end
    if size(ws.dv) != (nt, na)
        ws.dv       = Matrix{Float64}(undef, nt, na)
        ws.c3       = Matrix{Float64}(undef, nt, na)
        ws.vinf_arr = Matrix{Float64}(undef, nt, na)
        ws.ok       = BitMatrix(undef, nt, na)
    end
    return ws
end

@inline function dv_from_vinf(μp::Float64, rpark::Float64, vinf::Float64)
    vcirc = sqrt(μp / rpark)
    vesc  = sqrt(2*μp / rpark)
    return sqrt(vinf*vinf + vesc*vesc) - vcirc
end

function precompute_states!(rs::Vector{Vec3}, vs::Vector{Vec3},
                            sys::TwoBodySystem, bm::BodyModel, ts::Vector{Float64})
    @inbounds for k in eachindex(ts)
        rs[k], vs[k] = state(sys, bm, ts[k])
    end
    return rs, vs
end

function porkchop_grid_dep_arr!(ws::PorkchopWorkspace,
                                sys::TwoBodySystem, dep::BodyModel, arr::BodyModel,
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

    ensure_workspace!(ws, nt, na)

    @inbounds for i in 1:nt
        ws.tdep[i] = t0 + (i-1)*dt
    end
    @inbounds for j in 1:na
        ws.tarr[j] = a0 + (j-1)*da
    end

    μc = sys.central.μ
    use_parking = (metric === :dv) && isfinite(dep_rpark) && isfinite(arr_rpark)

    precompute_states!(ws.r1s, ws.v1ps, sys, dep, ws.tdep)
    precompute_states!(ws.r2s, ws.v2ps, sys, arr, ws.tarr)

    dv = ws.dv
    c3 = ws.c3
    va = ws.vinf_arr
    ok = ws.ok

    # KEY FIX 1: :static scheduling reduces per-call scheduling overhead
    Threads.@threads :static for i in 1:nt
        td  = ws.tdep[i]
        r1  = ws.r1s[i]
        v1p = ws.v1ps[i]

        @inbounds for j in 1:na
            ta  = ws.tarr[j]
            tof = ta - td

            if tof < tof_min || tof > tof_max
                ok[i,j] = false
                dv[i,j] = NaN
                c3[i,j] = NaN
                va[i,j] = NaN
                continue
            end

            r2  = ws.r2s[j]
            v2p = ws.v2ps[j]

            # KEY FIX 2: positional call in the hot loop (no keyword dispatch)
            v1t, v2t, success = lambert_uv(μc, r1, r2, tof, longway)
            ok[i,j] = success

            if !success
                dv[i,j] = NaN
                c3[i,j] = NaN
                va[i,j] = NaN
                continue
            end

            vdep_inf_vec = v1t - v1p
            varr_inf_vec = v2t - v2p

            c3ij     = vdep_inf_vec ⋅ vdep_inf_vec
            vinf_dep = sqrt(c3ij)
            vinfA    = sqrt(varr_inf_vec ⋅ varr_inf_vec)

            c3[i,j] = c3ij
            va[i,j] = vinfA

            dvij = if metric === :vinf
                vinf_dep + vinfA
            elseif use_parking
                dv_from_vinf(dep_μ, dep_rpark, vinf_dep) + dv_from_vinf(arr_μ, arr_rpark, vinfA)
            else
                vinf_dep + vinfA
            end

            if dvij > dv_cap
                ok[i,j] = false
                dv[i,j] = NaN
            else
                dv[i,j] = dvij
            end
        end
    end

    return PorkchopResult(ws.tdep, ws.tarr, dv, c3, va, ok)
end

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

    ws = PorkchopWorkspace(nt, na)
    return porkchop_grid_dep_arr!(ws, sys, dep, arr, tdep_range, tarr_range;
        longway=longway, metric=metric,
        dep_μ=dep_μ, dep_rpark=dep_rpark,
        arr_μ=arr_μ, arr_rpark=arr_rpark,
        tof_min=tof_min, tof_max=tof_max, dv_cap=dv_cap
    )
end