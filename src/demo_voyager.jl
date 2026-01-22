include("PorkchopSolver.jl")
using .PorkchopSolver
using Dates
using Downloads
using CairoMakie

# ==============================================================================
# 1. SETUP AMBIENTE & KERNELS
# ==============================================================================
# Percorso relativo: cartella "kernels" sorella di "src"
kern_dir = joinpath(@__DIR__, "..", "kernels")
mkpath(kern_dir)
PS = PorkchopSolver

println("Checking kernels in: $kern_dir")
kernels = Dict(
    "naif0012.tls" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
    "pck00010.tpc" => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc",
    "de440.bsp"    => "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
)

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

# ==============================================================================
# 2. DEFINIZIONE SISTEMA SOLARE
# ==============================================================================
μ_sun = 1.32712440018e11
sun = PS.Body("Sun", μ_sun, 695700.0)
sys = PS.TwoBodySystem(sun)

earth   = PS.Body("Earth", 3.986004354e5, 6378.1363)
jupiter = PS.Body("Jupiter", 1.26686534e8, 71492.0)
saturn  = PS.Body("Saturn", 3.7931187e7, 60268.0)

# Ephemeris Models
ephE = PS.SpiceEphemeris("EARTH BARYCENTER", "SUN", "J2000", "NONE")
ephJ = PS.SpiceEphemeris("JUPITER BARYCENTER", "SUN", "J2000", "NONE")
ephS = PS.SpiceEphemeris("SATURN BARYCENTER", "SUN", "J2000", "NONE")

bm_dep = PS.BodyModel(earth, ephE)
bm_fly = PS.BodyModel(jupiter, ephJ)
bm_arr = PS.BodyModel(saturn, ephS)

# ==============================================================================
# 3. FINESTRA DI LANCIO VOYAGER 2 (1977)
# ==============================================================================
dep_start = "1977-07-20T00:00:00" 
dep_end   = "1977-10-10T00:00:00"

fly_start = "1979-06-15T00:00:00"
fly_end   = "1979-08-01T00:00:00"

arr_start = "1981-08-01T00:00:00"
arr_end   = "1981-09-15T00:00:00"

t0_d = PS.utc_to_et(dep_start); t1_d = PS.utc_to_et(dep_end)
t0_f = PS.utc_to_et(fly_start); t1_f = PS.utc_to_et(fly_end)
t0_a = PS.utc_to_et(arr_start); t1_a = PS.utc_to_et(arr_end)

# Risoluzione alta per il grafico (1 giorno)
dt = 1.0 * 86400.0

println("\n=== VOYAGER 2 VALIDATION RUN ===")
println("Target: Earth -> Jupiter -> Saturn")

# ==============================================================================
# 4. ESECUZIONE SOLVER
# ==============================================================================
println("Running optimized porkchop solver...")

# Eseguiamo direttamente il solver. 
@time res = PS.porkchop_flyby(sys, bm_dep, bm_fly, bm_arr,
    (t0_d, t1_d, dt),
    (t0_f, t1_f, dt),
    (t0_a, t1_a, dt);
    R_fly_min = jupiter.radius + 100000.0,
    rpark_dep = 6678.0, 
    rpark_arr = NaN, # Flyby a Saturno
    tof1_min  = 500 * 86400.0, tof1_max  = 800 * 86400.0,
    tof2_min  = 500 * 86400.0, tof2_max  = 1000 * 86400.0,
    dv_cap    = Inf
)

# ==============================================================================
# 5. ANALISI RISULTATI
# ==============================================================================
min_dv = Inf
best_idx = (0,0)
for i in axes(res.dv_total, 1), k in axes(res.dv_total, 2)
    val = res.dv_total[i,k]
    if isfinite(val) && val < min_dv
        global min_dv = val
        global best_idx = (i,k)
    end
end

if min_dv == Inf
    println("[!] No solution found.")
else
    i_opt, k_opt = best_idx
    
    # Estraiamo i tempi ottimali per la visualizzazione 3D
    t_dep_opt = res.tdep[i_opt]
    t_fly_opt = res.best_tfly[i_opt, k_opt]
    t_arr_opt = res.tarr[k_opt]

    println("\n=== BEST TRAJECTORY FOUND ===")
    println("  Departure (Earth):   ", PS.et_to_utc(t_dep_opt; prec=0))
    println("  Flyby (Jupiter):     ", PS.et_to_utc(t_fly_opt; prec=0))
    println("  Arrival (Saturn):    ", PS.et_to_utc(t_arr_opt; prec=0))
    println("  Total Delta V:       $(round(min_dv, digits=4)) km/s")

    # ==========================================================================
    # 6. PLOT 1: PORKCHOP MAP (2D)
    # ==========================================================================
    println("\nGenerating Porkchop Plot...")
    epoch = Date(dep_start[1:10])

    floor_min = floor(min_dv * 10) / 10.0
    fine_levels = collect(floor_min : 0.2 : 20.0)
    coarse_levels = collect(21.0 : 1.0 : 25.0)
    my_levels = unique(vcat(fine_levels, coarse_levels))
    
    fig_pork = PS.plot_grid(res;
        t0=t0_d, epoch0=epoch,
        Z=res.dv_total,
        title="Voyager 2 Recreation (E-J-S) - High Detail",
        zlabel="Total dV (km/s)",
        contour_levels=my_levels,
        zrange=(floor_min, 25.0),
        colormap=Reverse(:plasma)
    )

    pork_path = abspath(joinpath(@__DIR__, "..", "voyager_porkchop.png"))
    save(pork_path, fig_pork)
    println("Porkchop plot saved to: $pork_path")

    # ==========================================================================
    # 7. PLOT 2: TRAIETTORIA 3D (NUOVO!)
    # ==========================================================================
    println("Generating 3D Trajectory Visualization...")
    
    fig_traj = PS.plot_trajectory_3d(
        sys, 
        (bm_dep, bm_fly, bm_arr), 
        (t_dep_opt, t_fly_opt, t_arr_opt);
        title_str = "Voyager 2: Earth -> Jupiter -> Saturn"
    )
    
    traj_path = abspath(joinpath(@__DIR__, "..", "voyager_trajectory_3d.png"))
    save(traj_path, fig_traj)
    println("3D Trajectory saved to: $traj_path")
end