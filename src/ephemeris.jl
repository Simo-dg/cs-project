# src/ephemeris.jl
# Ephemeris backends: circular, Kepler, (optional) NAIF SPICE
# - For SPICE, we keep calls thread-safe by using a lock here.
# - For performance, porkchop precomputes states serially anyway.

using StaticArrays
using LinearAlgebra
using Base: ReentrantLock

# -----------------------------------------------------------------------------
# Circular coplanar ephemeris (fast demo model)
# -----------------------------------------------------------------------------
@inline function state(eph::CircularCoplanarEphemeris, t::Float64)
    θ = eph.λ0 + eph.n*t
    c = cos(θ); s = sin(θ)
    r = Vec3(eph.a*c, eph.a*s, 0.0)
    v = Vec3(-eph.a*eph.n*s, eph.a*eph.n*c, 0.0)
    return r, v
end

# -----------------------------------------------------------------------------
# Kepler ephemeris (two-body around system.central)
# -----------------------------------------------------------------------------
function state(μcentral::Float64, eph::KeplerEphemeris, t::Float64)
    el = eph.el
    dt = t - el.t0
    n = sqrt(μcentral/(el.a^3))
    M = el.M0 + n*dt
    el2 = KeplerianElements(el.a, el.e, el.i, el.Ω, el.ω, M, el.t0)
    return kepler_to_rv(μcentral, el2)
end

# -----------------------------------------------------------------------------
# Generic: state(sys, BodyModel, t)
# -----------------------------------------------------------------------------
@inline function state(sys::TwoBodySystem, bm::BodyModel, t::Float64)
    if bm.eph isa CircularCoplanarEphemeris
        return state(bm.eph, t)
    elseif bm.eph isa KeplerEphemeris
        return state(sys.central.μ, bm.eph, t)
    else
        return state(sys, bm.eph, t)  # e.g., SPICE below
    end
end

# -----------------------------------------------------------------------------
# Optional SPICE backend (NAIF kernels, e.g. de440.bsp)
# -----------------------------------------------------------------------------
const _HAS_SPICE = let ok = false
    try
        @eval import SPICE
        ok = true
    catch
        ok = false
    end
    ok
end

const _SPICE_LOCK = ReentrantLock()
const _SSB = "SOLAR SYSTEM BARYCENTER"

struct SpiceEphemeris <: AbstractEphemeris
    target::String    # e.g. "EARTH BARYCENTER", "MARS BARYCENTER"
    observer::String  # use "SUN" for heliocentric (implemented via subtraction), or _SSB for barycentric
    frame::String     # "J2000"
    abcorr::String    # "NONE"
end

function spice_load_kernels!(paths::AbstractVector{<:AbstractString})
    _HAS_SPICE || error("SPICE.jl not installed. Run: julia --project -e 'import Pkg; Pkg.add(\"SPICE\")'")
    lock(_SPICE_LOCK)
    try
        for p in paths
            SPICE.furnsh(String(p))
        end
    finally
        unlock(_SPICE_LOCK)
    end
    return nothing
end

function utc_to_et(utc::AbstractString)
    _HAS_SPICE || error("SPICE.jl not installed. Run: julia --project -e 'import Pkg; Pkg.add(\"SPICE\")'")
    lock(_SPICE_LOCK)
    try
        return SPICE.utc2et(String(utc))
    finally
        unlock(_SPICE_LOCK)
    end
end

@inline function _spkezr_locked(target::String, et::Float64, frame::String, abcorr::String, obs::String)
    st, _lt = SPICE.spkezr(target, et, frame, abcorr, obs)
    return st
end

function state(::TwoBodySystem, eph::SpiceEphemeris, et::Float64)
    _HAS_SPICE || error("SPICE.jl is not installed, cannot use SpiceEphemeris.")
    lock(_SPICE_LOCK)
    try
        if eph.observer == "SUN"
            # Robust heliocentric: target wrt SSB minus Sun wrt SSB
            stT = _spkezr_locked(eph.target, et, eph.frame, eph.abcorr, _SSB)
            stS = _spkezr_locked("SUN",     et, eph.frame, eph.abcorr, _SSB)
            r = Vec3(stT[1]-stS[1], stT[2]-stS[2], stT[3]-stS[3])
            v = Vec3(stT[4]-stS[4], stT[5]-stS[5], stT[6]-stS[6])
            return r, v
        else
            st = _spkezr_locked(eph.target, et, eph.frame, eph.abcorr, eph.observer)
            r = Vec3(st[1], st[2], st[3])
            v = Vec3(st[4], st[5], st[6])
            return r, v
        end
    finally
        unlock(_SPICE_LOCK)
    end
end