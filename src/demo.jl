# src/demo.jl
# Run:
#   julia --project -t auto -e 'include("src/demo.jl")'
#
# Installs (once):
#   julia --project -e 'import Pkg; Pkg.add(["SPICE","Downloads","CairoMakie","StaticArrays"])'

include("PorkchopSolver.jl")
using .PorkchopSolver
using Dates
using CairoMakie
using Downloads

const PS = PorkchopSolver

function main()
    # -------------------------------------------------
    # 1) Ensure NAIF kernels exist (auto-download)
    # -------------------------------------------------
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
        if !isfile(fpath)
            println("  downloading $fname")
            Downloads.download(url, fpath)
        else
            println("  found $fname")
        end
    end

    PS.spice_load_kernels!([
        joinpath(kern_dir, "naif0012.tls"),
        joinpath(kern_dir, "pck00010.tpc"),
        joinpath(kern_dir, "de440.bsp"),
    ])

    # -------------------------------------------------
    # 2) Two-body Lambert central body (Sun)
    # -------------------------------------------------
    μ_sun = 1.32712440018e11  # km^3/s^2
    sun = Body("Sun", μ_sun, 695700.0)
    sys = TwoBodySystem(sun)

    earth = Body("Earth", 3.986004354e5, 6378.1363)
    mars  = Body("Mars",  4.282837e4,   3396.19)

    # Use barycenters for reliability with DE kernels
    ephE = PS.SpiceEphemeris("EARTH BARYCENTER", "SUN", "J2000", "NONE")
    ephM = PS.SpiceEphemeris("MARS BARYCENTER",  "SUN", "J2000", "NONE")

    dep = BodyModel(earth, ephE)
    arr = BodyModel(mars,  ephM)

    # -------------------------------------------------
    # 3) Define UTC date ranges, convert to ET
    # -------------------------------------------------
    dep_start_utc = "2026-01-01T00:00:00"
    dep_end_utc   = "2028-01-01T00:00:00"

    arr_start_utc = "2026-03-01T00:00:00"
    arr_end_utc   = "2029-01-01T00:00:00"

    tdep0 = PS.utc_to_et(dep_start_utc)
    tdep1 = PS.utc_to_et(dep_end_utc)
    tarr0 = PS.utc_to_et(arr_start_utc)
    tarr1 = PS.utc_to_et(arr_end_utc)

    dt = 1 * 86400.0
    tdep_range = (tdep0, tdep1, dt)
    tarr_range = (tarr0, tarr1, dt)

    epoch0 = Date(2026, 1, 1)  # axis labels

    # -------------------------------------------------
    # 4) Mission assumptions
    # -------------------------------------------------
    rpark_earth = 6678.0
    rpark_mars  = 4500.0

    tof_min =  90 * 86400.0
    tof_max = 350 * 86400.0
    dv_cap  = 30.0

    # -------------------------------------------------
    # 5) Compute porkchop
    # -------------------------------------------------
    println("Computing porkchop grid...")
    res = porkchop_grid_dep_arr(sys, dep, arr, tdep_range, tarr_range;
        metric=:dv,
        dep_rpark=rpark_earth,
        arr_rpark=rpark_mars,
        tof_min=tof_min,
        tof_max=tof_max,
        dv_cap=dv_cap
    )

    # -------------------------------------------------
    # 6) Stats + best finite solution (no NaNs)
    # -------------------------------------------------
    finite = filter(isfinite, vec(res.dv))
    dvmin = minimum(finite)
    dvmax = maximum(finite)
    @show dvmin dvmax

    best_dv = Inf
    best_i  = 0
    best_j  = 0
    @inbounds for i in axes(res.dv, 1), j in axes(res.dv, 2)
        v = res.dv[i, j]
        if isfinite(v) && v < best_dv
            best_dv = v
            best_i  = i
            best_j  = j
        end
    end

    best_dep_et = res.tdep[best_i]
    best_arr_et = res.tarr[best_j]

    best_dep_days = round(Int, (best_dep_et - tdep0) / 86400.0)
    best_arr_days = round(Int, (best_arr_et - tdep0) / 86400.0)

    best_dep_date = epoch0 + Day(best_dep_days)
    best_arr_date = epoch0 + Day(best_arr_days)
    best_tof_days = round((best_arr_et - best_dep_et) / 86400; digits=1)

    println("\n--- BEST EARTH → MARS (SPICE DE440, finite-only min Δv) ---")
    @show best_dep_date best_arr_date best_tof_days best_dv

    # -------------------------------------------------
    # 7) Plot + save
    # -------------------------------------------------
    fig = plot_porkchop(res;
        epoch0=epoch0,
        title="Earth → Mars porkchop (SPICE DE440) — Δv total (km/s)"
    )

    CairoMakie.save("porkchop_spice_dep_arr_dv.png", fig)
    display(fig)
end

main()