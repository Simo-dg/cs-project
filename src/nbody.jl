using LinearAlgebra
using StaticArrays

# Robust definition: V must be a vector, and elements must be subtypes of BodyModel
struct NBodySystem{V <: AbstractVector}
    bodies::V
    function NBodySystem(bodies::V) where {V <: AbstractVector}
        new{V}(bodies)
    end
end

# 1. THE PHYSICS
function compute_acceleration(sys::NBodySystem, r_sc::Vec3, t::Float64)
    acc = Vec3(0.0, 0.0, 0.0)
    
    # We construct a temporary TwoBodySystem wrapper using the first body (Sun).
    # This allows us to use the unified 'state(sys, bm, t)' interface from ephemeris.jl
    # which correctly handles dispatch for SPICE vs Circular vs Kepler.
    # Note: sys.bodies[1].body is usually the Sun.
    tsys = TwoBodySystem(sys.bodies[1].body)
    
    @inbounds for i in eachindex(sys.bodies)
        bm = sys.bodies[i]
        
        # CORRECT CALL: Pass the TwoBodySystem wrapper and the BodyModel.
        # This routes to line 90 in ephemeris.jl automatically.
        r_planet, _ = state(tsys, bm, t) 
        
        d_vec = r_planet - r_sc
        dist_sq = dot(d_vec, d_vec)
        dist = sqrt(dist_sq)
        
        # Avoid singularity (e.g. if starting exactly at Earth)
        if dist < 6500.0 
            continue 
        end
        
        μ = bm.body.μ
        acc += (μ / (dist_sq * dist)) * d_vec
    end
    
    return acc
end

# 2. THE INTEGRATOR (Runge-Kutta 4)
function rk4_step(sys::NBodySystem, r::Vec3, v::Vec3, t::Float64, dt::Float64)
    a1 = compute_acceleration(sys, r, t)
    v1 = v
    
    r2 = r + v1 * (0.5 * dt)
    a2 = compute_acceleration(sys, r2, t + 0.5 * dt)
    v2 = v + a1 * (0.5 * dt)
    
    r3 = r + v2 * (0.5 * dt)
    a3 = compute_acceleration(sys, r3, t + 0.5 * dt)
    v3 = v + a2 * (0.5 * dt)
    
    r4 = r + v3 * dt
    a4 = compute_acceleration(sys, r4, t + dt)
    v4 = v + a3 * dt
    
    r_new = r + (dt / 6.0) * (v1 + 2*v2 + 2*v3 + v4)
    v_new = v + (dt / 6.0) * (a1 + 2*a2 + 2*a3 + a4)
    
    return r_new, v_new
end

# 3. THE PROPAGATOR LOOP
function propagate_nbody(sys::NBodySystem, r0::Vec3, v0::Vec3, t0::Float64, tf::Float64; dt::Float64=3600.0)
    n_steps = Int(ceil((tf - t0) / dt))
    
    times = Vector{Float64}(undef, n_steps + 1)
    pos   = Vector{Vec3}(undef, n_steps + 1)
    vel   = Vector{Vec3}(undef, n_steps + 1)
    
    times[1] = t0
    pos[1]   = r0
    vel[1]   = v0
    
    t_curr = t0
    r_curr = r0
    v_curr = v0
    
    for i in 1:n_steps
        r_next, v_next = rk4_step(sys, r_curr, v_curr, t_curr, dt)
        
        t_curr += dt
        r_curr = r_next
        v_curr = v_next
        
        times[i+1] = t_curr
        pos[i+1]   = r_curr
        vel[i+1]   = v_curr
    end
    
    return times, pos, vel
end