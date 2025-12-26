using LinearAlgebra
using StaticArrays

@inline function R3(θ::Float64)
    c = cos(θ); s = sin(θ)
    return SMatrix{3,3,Float64,9}( c, s, 0,
                                 -s, c, 0,
                                  0, 0, 1)
end

@inline function R1(θ::Float64)
    c = cos(θ); s = sin(θ)
    return SMatrix{3,3,Float64,9}(1, 0, 0,
                                  0, c, s,
                                  0,-s, c)
end

@inline function perifocal_to_eci(i::Float64, Ω::Float64, ω::Float64)
    return R3(Ω) * (R1(i) * R3(ω))
end

@inline function kepler_E(M::Float64, e::Float64)
    E = M
    for _ in 1:12
        f  = E - e*sin(E) - M
        fp = 1 - e*cos(E)
        dE = -f/fp
        E += dE
        if abs(dE) < 1e-13
            break
        end
    end
    return E
end

@inline function kepler_to_rv(μ::Float64, el)
    a,e,i,Ω,ω,M0 = el.a,el.e,el.i,el.Ω,el.ω,el.M0
    E = kepler_E(M0, e)
    cE = cos(E); sE = sin(E)

    r_p = a*(1 - e*cE)
    x = a*(cE - e)
    y = a*sqrt(1 - e*e)*sE
    r_pf = Vec3(x, y, 0.0)

    n = sqrt(μ/(a^3))
    vx = -a*n*sE / r_p
    vy =  a*n*sqrt(1 - e*e)*cE / r_p
    v_pf = Vec3(vx, vy, 0.0)

    Q = perifocal_to_eci(i, Ω, ω)
    r = Q*r_pf
    v = Q*v_pf
    return r, v
end

function propagate_universal(μ::Float64, r0::Vec3, v0::Vec3, dt::Float64)
    r0n = sqrt(r0⋅r0)
    vr0 = (r0⋅v0)/r0n
    v02 = v0⋅v0
    α = 2/r0n - v02/μ

    χ = if abs(α) > 1e-10
        sqrt(μ)*abs(α)*dt
    else
        sqrt(μ)*dt/r0n
    end

    for _ in 1:30
        z = α*χ*χ
        C = stumpff_C(z)
        S = stumpff_S(z)

        F = (r0n*vr0/sqrt(μ))*χ*χ*C + (1 - α*r0n)*χ^3*S + r0n*χ - sqrt(μ)*dt
        dF = (r0n*vr0/sqrt(μ))*χ*(1 - z*S) + (1 - α*r0n)*χ^2*C + r0n

        dχ = -F/dF
        χ += dχ
        if abs(dχ) < 1e-11
            break
        end
    end

    z = α*χ*χ
    C = stumpff_C(z)
    S = stumpff_S(z)

    f = 1 - (χ*χ/r0n)*C
    g = dt - (1/sqrt(μ))*χ^3*S
    r = f*r0 + g*v0
    rn = sqrt(r⋅r)

    fdot = (sqrt(μ)/(rn*r0n))*(α*χ^3*S - χ)
    gdot = 1 - (χ*χ/rn)*C
    v = fdot*r0 + gdot*v0

    return r, v
end