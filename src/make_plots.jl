using CairoMakie
using JLD2
using Dates
using Printf
using Downloads

include("PorkchopSolver.jl")
using .PorkchopSolver
const PS = PorkchopSolver

# ------------------------------------------------------------------
# 1. KERNEL LOADING (REQUIRED FOR 3D TRAJECTORIES)
# ------------------------------------------------------------------
kern_dir = joinpath(@__DIR__, "..", "kernels")
mkpath(kern_dir)

kernels = Dict(
    "naif0012.tls" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "pck00010.tpc" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc",
    "de440.bsp"    => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
)

for (fname, url) in kernels
    fpath = joinpath(kern_dir, fname)
    if !isfile(fpath)
        Downloads.download(url, fpath)
    end
end

PS.spice_load_kernels!([
    joinpath(kern_dir, "naif0012.tls"),
    joinpath(kern_dir, "pck00010.tpc"),
    joinpath(kern_dir, "de440.bsp"),
])

# ------------------------------------------------------------------
# 2. PLOTTING LOGIC
# ------------------------------------------------------------------
data_dir = joinpath(@__DIR__, "..")
plot_dir = joinpath(@__DIR__, "..", "plots")
mkpath(plot_dir)

println("=== AUTOMATED PLOT GENERATION ===\n")

# --- MARS MISSION ---
path_mars = joinpath(data_dir, "mission_data.jld2")
if isfile(path_mars)
    println(">>> Processing Mars Data...")
    file = jldopen(path_mars, "r")
    res   = file["res"]
    t0    = file["tdep0"]
    close(file)

    epoch = Date(2026, 1, 1)
    
    save(joinpath(plot_dir, "mars_2026_dv.png"), 
         PS.plot_dv_total(res; t0=t0, epoch0=epoch, levels=collect(5:0.5:15), zrange=(5,15)))
         
    save(joinpath(plot_dir, "mars_2026_c3.png"), 
         PS.plot_c3(res; t0=t0, epoch0=epoch, levels=[10,20,30,40], zrange=(0,50)))
         
    println("    [OK] Mars plots generated.")
else
    println("    [SKIP] Mars data not found")
end

# --- FLYBY 2020 ---
path_fly = joinpath(data_dir, "flyby_data.jld2")
if isfile(path_fly)
    println("\n>>> Processing Flyby 2020 Data...")
    file = jldopen(path_fly, "r")
    res = file["res"]
    t0  = file["t0_d"]
    close(file)

    min_val = minimum(filter(isfinite, res.dv_total))
    plot_min = max(0.0, floor(min_val)-1)
    
    fig = PS.plot_grid(res; t0=t0, epoch0=Date(2019,6,1), Z=res.dv_total,
        title="Earth-Jup-Sat (2020)", zlabel="Total dV",
        contour_levels=collect(plot_min:0.5:plot_min+10), zrange=(plot_min, plot_min+10))
        
    save(joinpath(plot_dir, "flyby_2020_porkchop.png"), fig)
    println("    [OK] Flyby 2020 plot generated.")
else
    println("    [SKIP] Flyby data not found")
end

# --- VOYAGER 2 ---
path_voy = joinpath(data_dir, "voyager2_data.jld2")
if isfile(path_voy)
    println("\n>>> Processing Voyager 2 Data...")
    file = jldopen(path_voy, "r")
    res = file["res"]
    t0  = file["t0_d"]
    best_idx = file["best_idx"]
    sys = file["sys"]
    bds = (file["bm_dep"], file["bm_fly"], file["bm_arr"])
    close(file)

    Z_custom = res.dv_dep .+ res.dv_fly
    fig_2d = PS.plot_grid(res; t0=t0, epoch0=Date(1977,8,15), Z=Z_custom,
        title="Voyager 2 Cost", zlabel="dV (km/s)",
        contour_levels=collect(6.5:0.25:9.0), zrange=(6.5, 9.0))
    save(joinpath(plot_dir, "voyager2_porkchop.png"), fig_2d)

    i, k = best_idx
    times = (res.tdep[i], res.best_tfly[i,k], res.tarr[k])
    fig_3d = PS.plot_trajectory_3d(sys, bds, times; title_str="Voyager 2 Trajectory")
    save(joinpath(plot_dir, "voyager2_3d.png"), fig_3d)

    println("    [OK] Voyager plots generated.")
else
    println("    [SKIP] Voyager data not found")
end

println("\nAll tasks completed. Plots are in: $plot_dir")