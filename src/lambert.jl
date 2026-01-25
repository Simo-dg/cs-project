"""
    lambert.jl

    This module provides functions for solving Lambert's problem, which determines the orbit connecting two points
    in space in a given time under a central gravitational field.

    Main features:
    - lambert_uv: computes initial and final velocity vectors for a given transfer
    - Support for both short-way and long-way transfers
    - Uses Stumpff functions for universal variable formulation
"""


using LinearAlgebra
using StaticArrays


function lambert_uv(μ::Float64, r1::Vec3, r2::Vec3, dt::Float64, longway::Bool)
    r1n = norm(r1)
    r2n = norm(r2)

    cosΔν = (r1 ⋅ r2) / (r1n * r2n)
    cosΔν = clamp(cosΔν, -1.0, 1.0)
    Δν = acos(cosΔν)
    if longway
        Δν = 2π - Δν
    end

    sinΔν = sin(Δν)
    if abs(sinΔν) < 1e-14
        return (Vec3(NaN,NaN,NaN), Vec3(NaN,NaN,NaN), false)
    end

    A = sinΔν * sqrt((r1n * r2n) / (1 - cosΔν))
    if !isfinite(A) || abs(A) < 1e-15
        return (Vec3(NaN,NaN,NaN), Vec3(NaN,NaN,NaN), false)
    end

    z = 0.0
    tol = 1e-9
    converged = false

    for i in 1:60
        C = stumpff_C(z)
        S = stumpff_S(z)

        y = r1n + r2n + (A * (z * S - 1)) / sqrt(C)
        
        
        if y < 0
            z += 0.1
            continue
        end

        χ = sqrt(y / C)
        dtz = (χ^3 * S + A * sqrt(y)) / sqrt(μ)

        F = dtz - dt
        
        # Check for convergence
        if abs(F) < tol
            converged = true
            break
        end

        # Derivatives of Stumpff functions
        dC = abs(z) < 1e-8 ? (-1/24 + z/360 - z^2/13440) : (1 / (2z)) * (1 - z * S - 2C)
        dS = abs(z) < 1e-8 ? (-1/120 + z/2520 - z^2/120960) : (1 / (2z)) * (C - 3S)

        sqrtC = sqrt(C)
        dy = A * ((S + z * dS) * sqrtC - (z * S - 1) * dC / (2 * sqrtC)) / C
        dχ = 0.5 * ((dy * C - y * dC) / (C^2)) / (sqrt(y / C))

        ddt = ((3 * χ^2 * dχ) * S + χ^3 * dS + A * (0.5 / sqrt(y)) * dy) / sqrt(μ)

        dz = -F / ddt
        if !isfinite(dz)
            dz = sign(F) * 0.1
        end
        
        dz = clamp(dz, -2.0, 2.0)
        z += dz

        if abs(dz) < 1e-10
            converged = true 
            break
        end
    end

    if !converged
        return (Vec3(NaN,NaN,NaN), Vec3(NaN,NaN,NaN), false)
    end

    # Re-calculate final parameters
    C = stumpff_C(z)
    S = stumpff_S(z)
    y = r1n + r2n + (A * (z * S - 1)) / sqrt(C)
    
    if y <= 0
        return (Vec3(NaN,NaN,NaN), Vec3(NaN,NaN,NaN), false)
    end

    f    = 1 - y / r1n
    g    = A * sqrt(y / μ)
    gdot = 1 - y / r2n

    if abs(g) < 1e-14
        return (Vec3(NaN,NaN,NaN), Vec3(NaN,NaN,NaN), false)
    end

    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g
    
    return (v1, v2, true)
end

# default short-way
lambert_uv(μ::Float64, r1::Vec3, r2::Vec3, dt::Float64) =
    lambert_uv(μ, r1, r2, dt, false)


lambert_uv(μ::Float64, r1::Vec3, r2::Vec3, dt::Float64; longway::Bool=false) =
    lambert_uv(μ, r1, r2, dt, longway)