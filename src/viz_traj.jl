# src/viz_traj.jl
using CairoMakie
using StaticArrays
using LinearAlgebra

# Helper per convertire i nostri Vec3 (SVector) in Point3f di Makie
_to_p3(v::Vec3) = Point3f(v[1], v[2], v[3])

# Helper per propagare la traiettoria della sonda (Leg)
function _propagate_leg_points(μ_sun::Float64, r0::Vec3, v0::Vec3, tof::Float64; n_steps=300)
    points = Vector{Point3f}(undef, n_steps + 1)
    dt = tof / n_steps
    
    # Propagazione kepleriana punto per punto
    for i in 0:n_steps
        r_at_t, _ = propagate_universal(μ_sun, r0, v0, i*dt)
        points[i+1] = _to_p3(r_at_t)
    end
    return points
end

# Helper per generare un'orbita planetaria COMPLETA (360 gradi) per riferimento visivo
function _planet_orbit_trace(sys::TwoBodySystem, bm::BodyModel, t_center::Float64)
    # Stimiamo il periodo orbitale: T = 2π * sqrt(a^3 / μ)
    # Se è SPICE, non abbiamo 'a' facilmente, quindi usiamo un periodo fisso ampio
    # Terra ~1 anno, Giove ~12, Saturno ~29. Prendiamo 30 anni per sicurezza.
    years_30 = 30 * 365.25 * 86400.0
    
    t_start = t_center - years_30 / 2.0
    t_end   = t_center + years_30 / 2.0
    n_pts   = 500 # Risoluzione curva
    
    times = range(t_start, t_end, length=n_pts)
    return [_to_p3(state(sys, bm.eph, t)[1]) for t in times]
end

"""
    plot_trajectory_3d(sys, bodies, epochs; title="Mission Trajectory")

Disegna la traiettoria con vista "Top-Down" e orbite complete.
"""
function plot_trajectory_3d(sys::TwoBodySystem, 
                            bodies::Tuple{BodyModel, BodyModel, BodyModel}, 
                            epochs::Tuple{Float64, Float64, Float64};
                            title_str="Optimal Trajectory (E-J-S)")

    μ_sun = sys.central.μ
    bm_dep, bm_fly, bm_arr = bodies
    t_dep, t_fly, t_arr = epochs
    
    # 1. Stati esatti agli eventi
    r_dep, _ = state(sys, bm_dep.eph, t_dep)
    r_fly, _ = state(sys, bm_fly.eph, t_fly)
    r_arr, _ = state(sys, bm_arr.eph, t_arr)

    # 2. Calcolo Archi di Lambert (Sonda)
    tof1 = t_fly - t_dep
    tof2 = t_arr - t_fly
    
    v_sc1_dep, _, ok1 = lambert_uv(μ_sun, r_dep, r_fly, tof1, false)
    v_sc2_dep, _, ok2 = lambert_uv(μ_sun, r_fly, r_arr, tof2, false)
    
    if !ok1 || !ok2
        println("Warning: Lambert failed plotting trajectory.")
        return Figure()
    end

    # 3. Generazione Punti Sonda
    leg1_pts = _propagate_leg_points(μ_sun, r_dep, v_sc1_dep, tof1)
    leg2_pts = _propagate_leg_points(μ_sun, r_fly, v_sc2_dep, tof2)

    # 4. Generazione Orbite Pianeti (Complete)
    orb_dep = _planet_orbit_trace(sys, bm_dep, t_dep)
    orb_fly = _planet_orbit_trace(sys, bm_fly, t_fly)
    orb_arr = _planet_orbit_trace(sys, bm_arr, t_arr)

    # 5. Creazione Scena
    # Usiamo una risoluzione quadrata per mantenere l'aspetto
    fig = Figure(size=(1000, 1000), backgroundcolor=:black)
    
    # Axis3 con vista fissata dall'alto (elevation = pi/2)
    ax = Axis3(fig[1, 1], 
        title=title_str, titlesize=24, titlecolor=:white,
        aspect=:data, # FONDAMENTALE: mantiene le proporzioni reali
        xlabel="x (km)", ylabel="y (km)", zlabel="z (km)",
        backgroundcolor=:black,
        xgridcolor=(:white, 0.1), ygridcolor=(:white, 0.1), zgridcolor=(:white, 0.1),
        xtickcolor=:white, ytickcolor=:white, ztickcolor=:white,
        xticklabelcolor=:white, yticklabelcolor=:white, zticklabelcolor=:transparent, # Nascondiamo Z labels
        azimuth = -pi/2, # Rotazione
        elevation = 0.9*pi/2, # Quasi perpendicolare (vista mappa)
        viewmode = :fit
    )

    # Sole
    meshscatter!(ax, [Point3f(0,0,0)], markersize=6e7, color=:yellow, label="Sun")

    # --- PIANETI ---
    # Orbite (linee sottili tratteggiate o trasparenti)
    lines!(ax, orb_dep, color=(:dodgerblue, 0.3), linewidth=1)
    lines!(ax, orb_fly, color=(:orange, 0.3), linewidth=1)
    lines!(ax, orb_arr, color=(:khaki, 0.3), linewidth=1) # Saturno color sabbia

    # Posizioni (pallini)
    meshscatter!(ax, [_to_p3(r_dep)], markersize=4e7, color=:dodgerblue, label="Earth Dep")
    meshscatter!(ax, [_to_p3(r_fly)], markersize=9e7, color=:orange, label="Jupiter Flyby")
    meshscatter!(ax, [_to_p3(r_arr)], markersize=8e7, color=:khaki, label="Saturn Arr")

    # --- SONDA ---
    # Leg 1: Terra -> Giove
    lines!(ax, leg1_pts, color=:cyan, linewidth=3, label="Leg 1")
    # Leg 2: Giove -> Saturno
    lines!(ax, leg2_pts, color=:magenta, linewidth=3, label="Leg 2")

    # Legenda
    axislegend(ax, position=:rt, backgroundcolor=(:black, 0.5), labelcolor=:white)

    return fig
end