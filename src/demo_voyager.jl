include("PorkchopSolver.jl")
using .PorkchopSolver
using Dates
using Downloads
using BenchmarkTools
using Profile
using JLD2

const PS = PorkchopSolver

# ------------------------------------------------------------------
# 1. SETUP
# ------------------------------------------------------------------
kern_dir = joinpath(@__DIR__, "..", "kernels")
mkpath(kern_dir)
# (Kernel download logic omitted for brevity, identical to flyby script)
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

# Voyager 2 Window (1977)
t0_d = PS.utc_to_et("1977-07-20T00:00:00")
t1_d = PS.utc_to_et("1977-10-10T00:00:00")
t0_f = PS.utc_to_et("1979-06-15T00:00:00")
t1_f = PS.utc_to_et("1979-08-01T00:00:00")
t0_a = PS.utc_to_et("1981-08-01T00:00:00")
t1_a = PS.utc_to_et("1981-09-15T00:00:00")

dt = 1.0 * 86400.0 # Resolution 1 day

println("\n=== VOYAGER 2 SIMULATION ===")

# ------------------------------------------------------------------
# 2. BENCHMARKING
# ------------------------------------------------------------------
voyager_task = () -> PS.porkchop_flyby(sys, bm_dep, bm_fly, bm_arr,
    (t0_d, t1_d, dt), (t0_f, t1_f, dt), (t0_a, t1_a, dt);
    R_fly_min = jupiter.radius + 100000.0,
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

println("\n=== @allocated Check ===")
GC.gc()
mem = @allocated voyager_task()
println("Memory: ", round(mem/1e6, digits=2), " MB")

# ------------------------------------------------------------------
# 3. EXTRACT OPTIMAL
# ------------------------------------------------------------------
min_dv = Inf
best_idx = (0,0)
for i in axes(res.dv_total,1), k in axes(res.dv_total,2)
    v = res.dv_total[i,k]
    if isfinite(v) && v < min_dv
        global min_dv = v
        global best_idx = (i,k)
    end
end
println("Best dV: $min_dv km/s")

# ------------------------------------------------------------------
# 4. SAVE
# ------------------------------------------------------------------
outfile = joinpath(@__DIR__, "..", "voyager_data.jld2")
jldsave(outfile; res, min_dv, best_idx, t0_d, sys, bm_dep, bm_fly, bm_arr)
println("Saved to: $outfile")