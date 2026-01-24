"""
    types.jl

    This module defines fundamental types and data structures for astrodynamics and interplanetary mission analysis.
    It includes representations for celestial bodies, two-body systems, orbital elements, ephemeris models, and body models.

    Main features:
    - Vec3: 3D vector type for positions and velocities
    - Body: physical properties of a celestial body
    - TwoBodySystem: central body system for orbit propagation
    - KeplerianElements: classical orbital elements
    - AbstractEphemeris and concrete ephemeris types (circular, Keplerian)
    - BodyModel: combines a body with its ephemeris
"""


using StaticArrays

const Vec3 = SVector{3,Float64}

struct Body
    name::String
    μ::Float64
    radius::Float64
end

struct TwoBodySystem
    central::Body
end

struct KeplerianElements
    a::Float64
    e::Float64
    i::Float64
    Ω::Float64
    ω::Float64
    M0::Float64
    t0::Float64
end

abstract type AbstractEphemeris end

struct CircularCoplanarEphemeris <: AbstractEphemeris
    a::Float64
    n::Float64
    λ0::Float64
end

struct KeplerEphemeris <: AbstractEphemeris
    el::KeplerianElements
end

struct BodyModel{E<:AbstractEphemeris}
    body::Body
    eph::E
end

@inline norm2(v::Vec3) = v ⋅ v