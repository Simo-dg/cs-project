"""
    zoom_analysis.jl

    This script performs a focused (zoomed-in) analysis of porkchop plots and transfer opportunities for a specific mission scenario.
    It loads the required SPICE kernels, sets up the planetary bodies and time windows, and generates high-resolution plots for selected transfer windows.

    Main features:
    - Loads and configures mission scenario for detailed analysis
    - Computes porkchop plots with fine time resolution
    - Visualizes results using CairoMakie

"""


include("PorkchopSolver.jl")
using .PorkchopSolver
using Dates
using CairoMakie


PS = PorkchopSolver
kern_dir = joinpath(@__DIR__, "..", "kernels")
PS.spice_load_kernels!([
    joinpath(kern_dir, "naif0012.tls"),
    joinpath(kern_dir, "pck00010.tpc"),
    joinpath(kern_dir, "de440.bsp"),
])


sun = PS.Body("Sun", 1.32712440018e11, 695700.0)
sys = PS.TwoBodySystem(sun)
earth   = PS.Body("Earth", 3.986004354e5, 6378.1363)
jupiter = PS.Body("Jupiter", 1.26686534e8, 71492.0)
saturn  = PS.Body("Saturn", 3.7931187e7, 60268.0)

bm_dep = PS.BodyModel(earth,   PS.SpiceEphemeris("EARTH BARYCENTER", "SUN", "J2000", "NONE"))
bm_fly = PS.BodyModel(jupiter, PS.SpiceEphemeris("JUPITER BARYCENTER", "SUN", "J2000", "NONE"))
bm_arr = PS.BodyModel(saturn,  PS.SpiceEphemeris("SATURN BARYCENTER", "SUN", "J2000", "NONE"))


dep_start = "2020-03-01T00:00:00"
dep_end   = "2020-04-10T00:00:00"

fly_start = "2021-06-01T00:00:00"
fly_end   = "2022-01-01T00:00:00"
arr_start = "2029-01-01T00:00:00"
arr_end   = "2030-06-01T00:00:00"

t0_d = PS.utc_to_et(dep_start); t1_d = PS.utc_to_et(dep_end)
t0_f = PS.utc_to_et(fly_start); t1_f = PS.utc_to_et(fly_end)
t0_a = PS.utc_to_et(arr_start); t1_a = PS.utc_to_et(arr_end)


dt = 86400.0 

println("Running High-Res Zoom (March 2020)...")
res = PS.porkchop_flyby(sys, bm_dep, bm_fly, bm_arr,
    (t0_d, t1_d, dt), (t0_f, t1_f, dt*5.0), (t0_a, t1_a, dt*5.0); 
    R_fly_min = jupiter.radius + 100000.0,
    rpark_dep = 6678.0, rpark_arr = NaN,
    tof1_min  = 200*86400.0, tof1_max = 800*86400.0,
    tof2_min  = 2000*86400.0, tof2_max = 4000*86400.0,
    dv_cap    = Inf
)

#PLOT & SAVE
fig = PS.plot_grid(res;
    t0=t0_d, epoch0=Date(dep_start[1:10]),
    Z=res.dv_total,
    title="Zoom: Earth->Jup->Sat (Mar 2020)",
    zlabel="Total dV (km/s)",
    zrange=(11.0, 14.0), 
    contour_levels=collect(11.0:0.25:14.0)
)

out_path = joinpath(@__DIR__, "..","plots","zoom_2020_highres.png")
save(out_path, fig)
println("Saved: $out_path")
