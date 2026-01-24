"""
    demo_voyager2.jl

    This script demonstrates the setup and execution of a porkchop plot analysis for the Voyager 2 mission scenario using the PorkchopSolver module.
    It loads the required SPICE kernels, defines the planetary bodies and their ephemerides, sets up the simulation time windows for the historical Voyager 2 launch, and benchmarks the porkchop_flyby solver.
    The script finds the best mission delta-v (dV) and saves the results to a JLD2 file for further analysis or visualization.

    Main steps:
    1. Download and load required SPICE kernels if not present.
    2. Define the Sun, Earth, Jupiter, and Saturn with their physical and ephemeris properties.
    3. Set up the departure, flyby, and arrival time windows for the Voyager 2 mission (1977-1981).
    4. Run the porkchop_flyby solver and benchmark its performance.
    5. Find the best mission delta-v and save the results to disk.

"""


include("PorkchopSolver.jl")
using .PorkchopSolver
using Dates
using Downloads
using BenchmarkTools
using Profile
using JLD2

const PS = PorkchopSolver

# ------------------------------------------------------------------
# 1. SETUP & KERNELS
# ------------------------------------------------------------------
kern_dir = joinpath(@__DIR__, "..", "kernels")
mkpath(kern_dir)

kernels = Dict(
    "naif0012.tls" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "pck00010.tpc" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc",
    "de440.bsp"    => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
)

println("Checking Kernels...")
for (fname, url) in kernels
    fpath = joinpath(kern_dir, fname)
    if !isfile(fpath) || filesize(fpath) < 1000
        Downloads.download(url, fpath)
    end
end

PS.spice_load_kernels!([
    joinpath(kern_dir, "naif0012.tls"),
    joinpath(kern_dir, "pck00010.tpc"),
    joinpath(kern_dir, "de440.bsp"),
])

# ------------------------------------------------------------------
# 2. SYSTEM DEFINITION (CONST for Performance!)
# ------------------------------------------------------------------
const sun = PS.Body("Sun", 1.32712440018e11, 695700.0)
const sys = PS.TwoBodySystem(sun)

const earth   = PS.Body("Earth", 3.986004354e5, 6378.1363)
const jupiter = PS.Body("Jupiter", 1.26686534e8, 71492.0)
const saturn  = PS.Body("Saturn", 3.7931187e7, 60268.0)

const bm_dep = PS.BodyModel(earth,   PS.SpiceEphemeris("EARTH BARYCENTER", "SUN", "J2000", "NONE"))
const bm_fly = PS.BodyModel(jupiter, PS.SpiceEphemeris("JUPITER BARYCENTER", "SUN", "J2000", "NONE"))
const bm_arr = PS.BodyModel(saturn,  PS.SpiceEphemeris("SATURN BARYCENTER", "SUN", "J2000", "NONE"))

# ------------------------------------------------------------------
# 3. WINDOW CONFIGURATION (Sept 1977)
# ------------------------------------------------------------------
const t0_d = PS.utc_to_et("1977-08-15T00:00:00")
const t1_d = PS.utc_to_et("1977-09-15T00:00:00")

const t0_f = PS.utc_to_et("1979-06-01T00:00:00")
const t1_f = PS.utc_to_et("1979-08-01T00:00:00")

const t0_a = PS.utc_to_et("1981-08-01T00:00:00")
const t1_a = PS.utc_to_et("1981-10-01T00:00:00")

const dt = 1.0 * 86400.0 # 1 Day resolution

println("\n=== VOYAGER 2 SIMULATION (1977) ===")

# ------------------------------------------------------------------
# 4. BENCHMARK & EXECUTION
# ------------------------------------------------------------------
# This closure captures 'const' globals, so it is type-stable (Zero Allocations)
voyager_task = () -> PS.porkchop_flyby(sys, bm_dep, bm_fly, bm_arr,
    (t0_d, t1_d, dt), (t0_f, t1_f, dt), (t0_a, t1_a, dt);
    R_fly_min = jupiter.radius + 50000.0,
    rpark_dep = 6678.0, 
    rpark_arr = NaN, 
    tof1_min  = 500 * 86400.0, tof1_max  = 800 * 86400.0,
    tof2_min  = 500 * 86400.0, tof2_max  = 1000 * 86400.0,
    dv_cap    = Inf
)

println("\n=== WARMUP ===")
voyager_task() 

println("\n=== @time RUN ===")
@time res = voyager_task()

println("\n=== @allocated Check (Should be low!) ===")
GC.gc()
mem = @allocated voyager_task()
println("Memory: ", round(mem/1e6, digits=2), " MB")

# ------------------------------------------------------------------
# 5. FIND OPTIMAL MISSION & ANALYZE DATES
# ------------------------------------------------------------------
min_mission_dv = Inf
best_idx = (0,0)

for i in axes(res.dv_total, 1), k in axes(res.dv_total, 2)
    # Voyager 2 Cost: Launch + Flyby (Powered assist) only
    # Arrival at Saturn is often considered a flyby for Voyager, 
    # so we focus on the cost to get there and the assist maneuver.
    cost = res.dv_dep[i,k] + res.dv_fly[i,k]
    
    if isfinite(cost) && cost < min_mission_dv
        global min_mission_dv = cost
        global best_idx = (i,k)
    end
end

# Extract indices
i_best, k_best = best_idx

# Get optimal Epochs (ET)
dep_et = res.tdep[i_best] 
fly_et = res.best_tfly[i_best, k_best] 
arr_et = res.tarr[k_best]

# Convert to UTC Strings for readability
dep_utc = PS.et_to_utc(dep_et) 
fly_utc = PS.et_to_utc(fly_et) 
arr_utc = PS.et_to_utc(arr_et)

# Calculate Time of Flight (TOF) in Days
tof_earth_jup = (fly_et - dep_et) / 86400.0
tof_jup_sat   = (arr_et - fly_et) / 86400.0
tof_total     = (arr_et - dep_et) / 86400.0

println("\n--- OPTIMAL VOYAGER 2 TRAJECTORY (Calculated) ---")
println("Departure (Earth):  $dep_utc")
println("Flyby (Jupiter):    $fly_utc")
println("Arrival (Saturn):   $arr_utc")

println("\n--- MISSION DURATION ---")
println("Earth -> Jupiter:   $(round(tof_earth_jup, digits=2)) days")
println("Jupiter -> Saturn:  $(round(tof_jup_sat, digits=2)) days")
println("Total Duration:     $(round(tof_total/365.25, digits=2)) years")

println("\n--- DELTA-V BREAKDOWN (km/s) ---")
println("Launch dV (C3 proxy): $(round(res.dv_dep[i_best, k_best], digits=3))") 
println("Flyby dV (Maneuver):  $(round(res.dv_fly[i_best, k_best], digits=3))") 
println("Total Mission dV:     $(round(min_mission_dv, digits=3))")

# ------------------------------------------------------------------
# 6. SAVE DATA
# ------------------------------------------------------------------
outfile = joinpath(@__DIR__, "..", "voyager2_data.jld2")
jldsave(outfile; res, min_mission_dv, best_idx, t0_d, sys, bm_dep, bm_fly, bm_arr) 
println("\nData saved to: $outfile")