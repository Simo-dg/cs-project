# src/make_all_plots.jl
using CairoMakie
using JLD2
using Dates
using Printf

# Load the Solver Module (needed for Struct definitions)
include("PorkchopSolver.jl")
using .PorkchopSolver
const PS = PorkchopSolver

# Define where data lives
data_dir = joinpath(@__DIR__, "..")
plot_dir = joinpath(@__DIR__, "..", "plots")
mkpath(plot_dir)

println("=== AUTOMATED PLOT GENERATION ===\n")

# -------------------------------------------------------
# 1. MARS MISSION (from demo.jl -> mission_data.jld2)
# -------------------------------------------------------
path_mars = joinpath(data_dir, "mission_data.jld2")
if isfile(path_mars)
    println(">>> Processing Mars Data...")
    file = jldopen(path_mars, "r")
    res   = file["res"]
    t0    = file["tdep0"]
    close(file)

    epoch = Date(2026, 1, 1)
    
    # Save standard triplet
    save(joinpath(plot_dir, "mars_2026_dv.png"), 
         PS.plot_dv_total(res; t0=t0, epoch0=epoch, levels=collect(5:0.5:15), zrange=(5,15)))
         
    save(joinpath(plot_dir, "mars_2026_c3.png"), 
         PS.plot_c3(res; t0=t0, epoch0=epoch, levels=[10,20,30,40], zrange=(0,50)))
         
    println("    [OK] Mars plots generated.")
else
    println("    [SKIP] Mars data not found (run demo.jl)")
end

# -------------------------------------------------------
# 2. FLYBY 2020 (from demo_flyby.jl -> flyby_data.jld2)
# -------------------------------------------------------
path_fly = joinpath(data_dir, "flyby_data.jld2")
if isfile(path_fly)
    println("\n>>> Processing Flyby 2020 Data...")
    file = jldopen(path_fly, "r")
    res = file["res"]
    t0  = file["t0_d"]
    close(file)

    # Focus on efficient valley
    min_val = minimum(filter(isfinite, res.dv_total))
    plot_min = max(0.0, floor(min_val)-1)
    
    fig = PS.plot_grid(res; t0=t0, epoch0=Date(2019,6,1), Z=res.dv_total,
        title="Earth-Jup-Sat (2020)", zlabel="Total dV",
        contour_levels=collect(plot_min:0.5:plot_min+10), zrange=(plot_min, plot_min+10))
        
    save(joinpath(plot_dir, "flyby_2020_porkchop.png"), fig)
    println("    [OK] Flyby 2020 plot generated.")
else
    println("    [SKIP] Flyby data not found (run demo_flyby.jl)")
end

# -------------------------------------------------------
# 3. VOYAGER 2 (from demo_voyager2.jl -> voyager2_data.jld2)
# -------------------------------------------------------
path_voy = joinpath(data_dir, "voyager2_data.jld2")
if isfile(path_voy)
    println("\n>>> Processing Voyager 2 Data...")
    file = jldopen(path_voy, "r")
    res = file["res"]
    t0  = file["t0_d"]
    best_idx = file["best_idx"]
    # Load bodies/sys for 3D plot
    sys = file["sys"]
    bds = (file["bm_dep"], file["bm_fly"], file["bm_arr"])
    close(file)

    # 2D Plot (Custom Cost)
    Z_custom = res.dv_dep .+ res.dv_fly
    fig_2d = PS.plot_grid(res; t0=t0, epoch0=Date(1977,8,15), Z=Z_custom,
        title="Voyager 2 Cost", zlabel="dV (km/s)",
        contour_levels=collect(6.5:0.25:9.0), zrange=(6.5, 9.0))
    save(joinpath(plot_dir, "voyager2_porkchop.png"), fig_2d)

    # 3D Plot
    i, k = best_idx
    times = (res.tdep[i], res.best_tfly[i,k], res.tarr[k])
    fig_3d = PS.plot_trajectory_3d(sys, bds, times; title_str="Voyager 2 Trajectory")
    save(joinpath(plot_dir, "voyager2_3d.png"), fig_3d)

    println("    [OK] Voyager plots generated.")
else
    println("    [SKIP] Voyager data not found (run demo_voyager2.jl)")
end

println("\nAll tasks completed. Plots are in: $plot_dir")