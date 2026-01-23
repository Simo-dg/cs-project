include("PorkchopSolver.jl")
using .PorkchopSolver
using Dates
using Downloads
using BenchmarkTools
using Profile
using JLD2

const PS = PorkchopSolver

# ------------------------------------------------------------------
# 1. ENVIRONMENT & KERNELS
# ------------------------------------------------------------------
kern_dir = joinpath(@__DIR__, "..", "kernels")
mkpath(kern_dir)

kernels = Dict(
    "naif0012.tls" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "pck00010.tpc" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc",
    "de440.bsp"    => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
)

println("Checking NAIF kernels...")
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
# 2. SETUP SIMULATION
# ------------------------------------------------------------------
sun = PS.Body("Sun", 1.32712440018e11, 695700.0)
sys = PS.TwoBodySystem(sun)

earth   = PS.Body("Earth", 3.986004354e5, 6378.1363)
jupiter = PS.Body("Jupiter", 1.26686534e8, 71492.0)
saturn  = PS.Body("Saturn", 3.7931187e7, 60268.0)

# SPICE Ephemerides
bm_dep = PS.BodyModel(earth,   PS.SpiceEphemeris("EARTH BARYCENTER", "SUN", "J2000", "NONE"))
bm_fly = PS.BodyModel(jupiter, PS.SpiceEphemeris("JUPITER BARYCENTER", "SUN", "J2000", "NONE"))
bm_arr = PS.BodyModel(saturn,  PS.SpiceEphemeris("SATURN BARYCENTER", "SUN", "J2000", "NONE"))

# 2020 Window (Great Conjunction alignment)
t0_d = PS.utc_to_et("2019-06-01T00:00:00")
t1_d = PS.utc_to_et("2021-06-01T00:00:00")

t0_f = PS.utc_to_et("2021-01-01T00:00:00")
t1_f = PS.utc_to_et("2024-01-01T00:00:00")

t0_a = PS.utc_to_et("2026-01-01T00:00:00")
t1_a = PS.utc_to_et("2030-01-01T00:00:00")

dt = 5.0 * 86400.0 

println("\n=== SIMULATION SETTINGS (Flyby 2020) ===")
println("Scan Resolution: 5.0 days")

# ------------------------------------------------------------------
# 3. BENCHMARKING WRAPPER
# ------------------------------------------------------------------

solver_task = () -> PS.porkchop_flyby(sys, bm_dep, bm_fly, bm_arr,
    (t0_d, t1_d, dt),
    (t0_f, t1_f, dt),
    (t0_a, t1_a, dt);
    R_fly_min = jupiter.radius + 200000.0,
    rpark_dep = 6678.0, 
    rpark_arr = NaN,
    tof1_min  = 200 * 86400.0, tof1_max  = 1000 * 86400.0,
    tof2_min  = 1000 * 86400.0, tof2_max  = 4000 * 86400.0,
    dv_cap    = Inf
)

# ------------------------------------------------------------------
# 4. EXECUTION & PROFILING
# ------------------------------------------------------------------
println("\n=== WARMUP ===")
solver_task() 

println("\n=== @time RUN ===")
@time res = solver_task()

println("\n=== @allocated (Zero-Alloc Check) ===")
GC.gc()
mem = @allocated solver_task()
println("Memory allocated: ", round(mem/1e6, digits=2), " MB")

# Optional: Full Benchmark (takes longer)
# println("\n=== BenchmarkTools ===")
# @btime ($solver_task)()

# ------------------------------------------------------------------
# 5. SAVE RESULTS
# ------------------------------------------------------------------
outfile = joinpath(@__DIR__, "..", "flyby_data.jld2")
println("\n=== SAVING DATA ===")
jldsave(outfile; res, t0_d, sys, bm_dep, bm_fly, bm_arr)
println("Saved to: $outfile")