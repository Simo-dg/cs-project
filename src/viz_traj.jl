using CairoMakie
using StaticArrays
using LinearAlgebra


_to_p3(v::Vec3) = Point3f(v[1], v[2], v[3])


function _propagate_leg_points(μ_sun::Float64, r0::Vec3, v0::Vec3, tof::Float64; n_steps=300)
    points = Vector{Point3f}(undef, n_steps + 1)
    dt = tof / n_steps
    
    for i in 0:n_steps
        r_at_t, _ = propagate_universal(μ_sun, r0, v0, i*dt)
        points[i+1] = _to_p3(r_at_t)
    end
    return points
end


function _planet_orbit_trace(sys::TwoBodySystem, bm::BodyModel, t_center::Float64)
    years_30 = 30 * 365.25 * 86400.0
    
    t_start = t_center - years_30 / 2.0
    t_end   = t_center + years_30 / 2.0
    n_pts   = 500 
    
    times = range(t_start, t_end, length=n_pts)
    return [_to_p3(state(sys, bm.eph, t)[1]) for t in times]
end


function plot_trajectory_3d(sys::TwoBodySystem, 
                            bodies::Tuple{BodyModel, BodyModel, BodyModel}, 
                            epochs::Tuple{Float64, Float64, Float64};
                            title_str="Optimal Trajectory (E-J-S)")

    μ_sun = sys.central.μ
    bm_dep, bm_fly, bm_arr = bodies
    t_dep, t_fly, t_arr = epochs
    
    r_dep, _ = state(sys, bm_dep.eph, t_dep)
    r_fly, _ = state(sys, bm_fly.eph, t_fly)
    r_arr, _ = state(sys, bm_arr.eph, t_arr)


    tof1 = t_fly - t_dep
    tof2 = t_arr - t_fly
    
    v_sc1_dep, _, ok1 = lambert_uv(μ_sun, r_dep, r_fly, tof1, false)
    v_sc2_dep, _, ok2 = lambert_uv(μ_sun, r_fly, r_arr, tof2, false)
    
    if !ok1 || !ok2
        println("Warning: Lambert failed plotting trajectory.")
        return Figure()
    end


    leg1_pts = _propagate_leg_points(μ_sun, r_dep, v_sc1_dep, tof1)
    leg2_pts = _propagate_leg_points(μ_sun, r_fly, v_sc2_dep, tof2)


    orb_dep = _planet_orbit_trace(sys, bm_dep, t_dep)
    orb_fly = _planet_orbit_trace(sys, bm_fly, t_fly)
    orb_arr = _planet_orbit_trace(sys, bm_arr, t_arr)


    fig = Figure(size=(1000, 1000), backgroundcolor=:black)
    
    # Axis3 
    ax = Axis3(fig[1, 1], 
        title=title_str, titlesize=24, titlecolor=:white,
        aspect=:data, 
        xlabel="x (km)", ylabel="y (km)", zlabel="z (km)",
        backgroundcolor=:black,
        xgridcolor=(:white, 0.1), ygridcolor=(:white, 0.1), zgridcolor=(:white, 0.1),
        xtickcolor=:white, ytickcolor=:white, ztickcolor=:white,
        xticklabelcolor=:white, yticklabelcolor=:white, zticklabelcolor=:transparent, 
        azimuth = -pi/2, 
        elevation = 0.9*pi/2,
        viewmode = :fit
    )


    meshscatter!(ax, [Point3f(0,0,0)], markersize=6e7, color=:yellow, label="Sun")


    lines!(ax, orb_dep, color=(:dodgerblue, 0.3), linewidth=1)
    lines!(ax, orb_fly, color=(:orange, 0.3), linewidth=1)
    lines!(ax, orb_arr, color=(:khaki, 0.3), linewidth=1) 


    meshscatter!(ax, [_to_p3(r_dep)], markersize=4e7, color=:dodgerblue, label="Earth Dep")
    meshscatter!(ax, [_to_p3(r_fly)], markersize=9e7, color=:orange, label="Jupiter Flyby")
    meshscatter!(ax, [_to_p3(r_arr)], markersize=8e7, color=:khaki, label="Saturn Arr")


    lines!(ax, leg1_pts, color=:cyan, linewidth=3, label="Leg 1")

    lines!(ax, leg2_pts, color=:magenta, linewidth=3, label="Leg 2")


    axislegend(ax, position=:rt, backgroundcolor=(:black, 0.5), labelcolor=:white)

    return fig
end