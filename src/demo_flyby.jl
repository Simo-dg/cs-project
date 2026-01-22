include("PorkchopSolver.jl")
using .PorkchopSolver
using Dates
using Downloads
using CairoMakie


# ------------------------------------------------------------------
# 1. ENVIRONMENT & KERNELS
# ------------------------------------------------------------------
kern_dir = joinpath(@__DIR__, "..", "kernels")
mkpath(kern_dir)
PS = PorkchopSolver

kernels = Dict(
    "naif0012.tls" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "pck00010.tpc" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc",
    "de440.bsp"    => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
)

println("Checking NAIF kernels in: $kern_dir")
for (fname, url) in kernels
    fpath = joinpath(kern_dir, fname)
    if !isfile(fpath) || filesize(fpath) < 1000
        println("  Downloading $fname...")
        Downloads.download(url, fpath)
    end
end

PS.spice_load_kernels!([
    joinpath(kern_dir, "naif0012.tls"),
    joinpath(kern_dir, "pck00010.tpc"),
    joinpath(kern_dir, "de440.bsp"),
])

# ------------------------------------------------------------------
# 2. PHYSICAL SYSTEM
# ------------------------------------------------------------------
μ_sun = 1.32712440018e11
sun = PS.Body("Sun", μ_sun, 695700.0)
sys = PS.TwoBodySystem(sun)

earth   = PS.Body("Earth", 3.986004354e5, 6378.1363)
jupiter = PS.Body("Jupiter", 1.26686534e8, 71492.0)
saturn  = PS.Body("Saturn", 3.7931187e7, 60268.0)

ephE = PS.SpiceEphemeris("EARTH BARYCENTER", "SUN", "J2000", "NONE")
ephJ = PS.SpiceEphemeris("JUPITER BARYCENTER", "SUN", "J2000", "NONE")
ephS = PS.SpiceEphemeris("SATURN BARYCENTER", "SUN", "J2000", "NONE")

bm_dep = PS.BodyModel(earth, ephE)
bm_fly = PS.BodyModel(jupiter, ephJ)
bm_arr = PS.BodyModel(saturn, ephS)

# ------------------------------------------------------------------
# 3. MISSION WINDOW (2020 DEMO CONFIGURATION)
# ------------------------------------------------------------------
# Ideally suited for the Great Conjunction alignment

# Departure: Late 2019 to 2021
dep_start = "2019-06-01T00:00:00"
dep_end   = "2021-06-01T00:00:00"

# Flyby: 2021 to 2024
fly_start = "2021-01-01T00:00:00"
fly_end   = "2024-01-01T00:00:00"

# Arrival: 2026 to 2030 (Fast transfer via gravity assist)
arr_start = "2026-01-01T00:00:00"
arr_end   = "2030-01-01T00:00:00"

t0_d = PS.utc_to_et(dep_start); t1_d = PS.utc_to_et(dep_end)
t0_f = PS.utc_to_et(fly_start); t1_f = PS.utc_to_et(fly_end)
t0_a = PS.utc_to_et(arr_start); t1_a = PS.utc_to_et(arr_end)

dt_days = 5.0
dt = dt_days * 86400.0

println("\n=== STARTING MULTI-GRAVITY ASSIST SCAN (DEMO) ===")
println("  Strategy:   Earth -> Jupiter -> Saturn (2020 Window)")
println("  Resolution: $dt_days days")

# ------------------------------------------------------------------
# 4. SOLVER EXECUTION
# ------------------------------------------------------------------
println("  Running High-Performance Solver...")

@time res = PS.porkchop_flyby(sys, bm_dep, bm_fly, bm_arr,
    (t0_d, t1_d, dt),
    (t0_f, t1_f, dt),
    (t0_a, t1_a, dt);
    R_fly_min = jupiter.radius + 200000.0,
    rpark_dep = 6678.0, 
    rpark_arr = NaN,
    tof1_min  = 200 * 86400.0,
    tof1_max  = 1000 * 86400.0,
    tof2_min  = 1000 * 86400.0,
    tof2_max  = 4000 * 86400.0,
    dv_cap    = Inf
)

# ------------------------------------------------------------------
# 5. ANALYSIS
# ------------------------------------------------------------------
min_dv = Inf
best_idx = (0,0)

for i in axes(res.dv_total, 1), k in axes(res.dv_total, 2)
    val = res.dv_total[i,k]
    if isfinite(val) && val < min_dv
        min_dv = val
        best_idx = (i,k)
    end
end

if min_dv == Inf
    println("\n[!] No solution found.")
else
    i, k = best_idx
    println("\n=== OPTIMAL TRAJECTORY FOUND ===")
    println("  Departure (Earth):   ", PS.et_to_utc(res.tdep[i]; prec=0))
    println("  Flyby (Jupiter):     ", PS.et_to_utc(res.best_tfly[i,k]; prec=0))
    println("  Arrival (Saturn):    ", PS.et_to_utc(res.tarr[k]; prec=0))
    println("  ---------------------------")
    println("  Total Delta V:       $(round(min_dv, digits=4)) km/s")
    println("    - Launch C3:       $(round(res.c3[i,k], digits=2)) km²/s²")
    println("    - Launch dV:       $(round(res.dv_dep[i,k], digits=4)) km/s")
    println("    - Flyby dV:        $(round(res.dv_fly[i,k]*1000, digits=1)) m/s")
end

# ------------------------------------------------------------------
# 6. PLOTTING & SAVING (OUTSIDE SRC)
# ------------------------------------------------------------------
if min_dv != Inf
    println("\nGenerating Plot...")
    
    # 2020 window is efficient, so we can focus the plot scale
    # on lower values (e.g., 6.0 to 12.0 km/s)
    floor_min = floor(min_dv)
    plot_min = max(0.0, floor_min - 1.0)
    plot_max = plot_min + 10.0 
    
    my_levels = collect(plot_min:0.5:plot_max)

    epoch = Date(dep_start[1:10])
    
    # 1. Capture Figure
    fig = PS.plot_grid(res;
        t0=t0_d, epoch0=epoch,
        Z=res.dv_total,
        title="Earth -> Jupiter -> Saturn (2020 Demo)",
        zlabel="Total dV (km/s)",
        contour_levels=my_levels,
        zrange=(plot_min, plot_max)
    )

    # 2. Save OUTSIDE src folder
    out_path = joinpath(@__DIR__, "..", "plots", "demo_2020_porkchop.png")
    
    save(out_path, fig)
    println("Plot saved to: $(abspath(out_path))")
end
