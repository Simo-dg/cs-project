include("PorkchopSolver.jl")
using .PorkchopSolver
using Dates
using Downloads
using BenchmarkTools
using Profile
using JLD2

# Toggle this to false to skip the Profile steps
const ENABLE_PROFILING = false

const PS = PorkchopSolver

function merge_min(a::PS.PorkchopResult, b::PS.PorkchopResult)
    dv  = similar(a.dv)
    c3  = similar(a.c3)
    va  = similar(a.vinf_arr)
    ok  = falses(size(a.ok)...)

    @inbounds for i in axes(dv,1), j in axes(dv,2)
        da = a.dv[i,j]
        db = b.dv[i,j]
        if isfinite(da) && (!isfinite(db) || da <= db)
            dv[i,j] = da;  c3[i,j] = a.c3[i,j];  va[i,j] = a.vinf_arr[i,j];  ok[i,j] = a.ok[i,j]
        elseif isfinite(db)
            dv[i,j] = db;  c3[i,j] = b.c3[i,j];  va[i,j] = b.vinf_arr[i,j];  ok[i,j] = b.ok[i,j]
        else
            dv[i,j] = NaN; c3[i,j] = NaN;       va[i,j] = NaN;               ok[i,j] = false
        end
    end
    return PS.PorkchopResult(a.tdep, a.tarr, dv, c3, va, ok)
end

# -------------------------------------------------
# 0) Precision knob: grid resolution 
# -------------------------------------------------
const dt_hours = 3.0            
const dt = dt_hours * 3600.0

# Detect allocation tracking 
opts = Base.JLOptions()
tracking_alloc = (getfield(opts, :tracked_path) != C_NULL) || (getfield(opts, :malloc_log) != 0)

# -------------------------------------------------
# 1) Ensure NAIF kernels exist (auto-download)
# -------------------------------------------------
const kern_dir = joinpath(@__DIR__, "..", "kernels")
mkpath(kern_dir)

const kernels = Dict(
    "naif0012.tls" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "pck00010.tpc" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc",
    "de440.bsp"    => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
)

println("Checking NAIF kernels...")
for (fname, url) in kernels
    fpath = joinpath(kern_dir, fname)
    if !isfile(fpath)
        println("  downloading $fname")
        Downloads.download(url, fpath)
    end
end

PS.spice_load_kernels!([
    joinpath(kern_dir, "naif0012.tls"),
    joinpath(kern_dir, "pck00010.tpc"),
    joinpath(kern_dir, "de440.bsp"),
])

# -------------------------------------------------
# 2) Two-body system (Sun)
# -------------------------------------------------
const μ_sun = 1.32712440018e11
const sun = PS.Body("Sun", μ_sun, 695700.0)
const sys = PS.TwoBodySystem(sun)

const earth = PS.Body("Earth", 3.986004354e5, 6378.1363)
const mars  = PS.Body("Mars",  4.282837e4,   3396.19)

const ephE = PS.SpiceEphemeris("EARTH BARYCENTER", "SUN", "J2000", "NONE")
const ephM = PS.SpiceEphemeris("MARS BARYCENTER",  "SUN", "J2000", "NONE")
const dep = PS.BodyModel(earth, ephE)
const arr = PS.BodyModel(mars,  ephM)

# -------------------------------------------------
# 3) Date ranges
# -------------------------------------------------
const dep_start_utc = "2026-01-01T00:00:00"
const dep_end_utc   = "2028-01-01T00:00:00"
const arr_start_utc = "2026-03-01T00:00:00"
const arr_end_utc   = "2029-01-01T00:00:00"

const tdep0 = PS.utc_to_et(dep_start_utc)
const tdep1 = PS.utc_to_et(dep_end_utc)
const tarr0 = PS.utc_to_et(arr_start_utc)
const tarr1 = PS.utc_to_et(arr_end_utc)

const tdep_range = (tdep0, tdep1, dt)
const tarr_range = (tarr0, tarr1, dt)

# -------------------------------------------------
# 4) Mission assumptions
# -------------------------------------------------
const rpark_earth = 6678.0
const rpark_mars  = 4500.0

const tof_min =  60 * 86400.0
const tof_max = 500 * 86400.0
const dv_cap  = 30.0

# -------------------------------------------------
# 5) Print grid size + expected work
# -------------------------------------------------
const nt = Int(floor((tdep1 - tdep0)/dt)) + 1
const na = Int(floor((tarr1 - tarr0)/dt)) + 1

# Optional stats (doesn't affect performance)
combos_raw = nt * na
band_width = Int(clamp(floor((tof_max - tof_min)/dt), 0, na))
combos_lambert_est = nt * band_width

println("\n=== GRID SETTINGS ===")
println("dt = $(dt_hours) hours")
println("Grid: $nt x $na")
println("Est. Solves: $(round(2*combos_lambert_est/1e6, digits=2)) million")

# -------------------------------------------------
# 6) Workspaces (MUST BE CONST)
# -------------------------------------------------
const wsS = PS.PorkchopWorkspace(nt, na)
const wsL = PS.PorkchopWorkspace(nt, na)

shortway = () -> PS.porkchop_grid_dep_arr!(wsS, sys, dep, arr, tdep_range, tarr_range;
    metric=:dv, dep_rpark=rpark_earth, arr_rpark=rpark_mars,
    tof_min=tof_min, tof_max=tof_max, dv_cap=dv_cap, longway=false
)

longway = () -> PS.porkchop_grid_dep_arr!(wsL, sys, dep, arr, tdep_range, tarr_range;
    metric=:dv, dep_rpark=rpark_earth, arr_rpark=rpark_mars,
    tof_min=tof_min, tof_max=tof_max, dv_cap=dv_cap, longway=true
)

# -------------------------------------------------
# 7) Warmup
# -------------------------------------------------
println("\n=== WARMUP (compile) ===")
shortway(); longway()

# -------------------------------------------------
# 8) Timing
# -------------------------------------------------
println("\n=== @time short-way ===")
@time resS = shortway()

println("\n=== @time long-way ===")
@time resL = longway()


if ENABLE_PROFILING && !tracking_alloc
    println("\n=== BenchmarkTools (@btime) short-way ===")
    @btime ($shortway)()
else
    println("\n=== Skipped Profiling (ENABLE_PROFILING = false) ===")
end

# -------------------------------------------------
# 9) Merge + best
# -------------------------------------------------

function find_best_mission(dv_matrix)
    best_val = Inf
    best_idx = (1, 1)

    @inbounds for i in CartesianIndices(dv_matrix)
        val = dv_matrix[i]
        if val < best_val
            best_val = val
            best_idx = Tuple(i)
        end
    end
    return best_val, best_idx
end

res = merge_min(resS, resL)

# Fast search using the function
best_dv, (best_i, best_j) = find_best_mission(res.dv)

# Extract times
best_dep_et = res.tdep[best_i]
best_arr_et = res.tarr[best_j]

best_dep_utc = PS.et_to_utc(best_dep_et; prec=0)
best_arr_utc = PS.et_to_utc(best_arr_et; prec=0)

best_tof_hours = round((best_arr_et - best_dep_et) / 3600; digits=2)
best_tof_days  = round(best_tof_hours / 24; digits=2)

println("\n--- BEST EARTH → MARS ---")
println("best_dep_utc   = $best_dep_utc")
println("best_arr_utc   = $best_arr_utc")
println("best_tof_days  = $best_tof_days  (hours = $best_tof_hours)")
println("best_dv        = $best_dv")

println("\n=== SAVING RESULTS ===")
output_file = joinpath(@__DIR__, "..", "mission_data.jld2")

jldsave(output_file; 
    res,             
    tdep0,          
    best_dv,
    best_dep_utc,
    best_arr_utc
)