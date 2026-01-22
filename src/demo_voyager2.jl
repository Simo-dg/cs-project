include("PorkchopSolver.jl")
using .PorkchopSolver
using Dates
using Downloads
using CairoMakie


# ------------------------------------------------------------------
# 1. SETUP
# ------------------------------------------------------------------
PS = PorkchopSolver
kern_dir = joinpath(@__DIR__, "..", "kernels")
mkpath(kern_dir)

# Caricamento Kernels
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
# 2. SISTEMA E CORPI
# ------------------------------------------------------------------
sun = PS.Body("Sun", 1.32712440018e11, 695700.0)
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
# 3. FINESTRA STORICA VOYAGER 2 (SETTEMBRE 1977)
# ------------------------------------------------------------------
dep_start = "1977-08-15T00:00:00"
dep_end   = "1977-09-15T00:00:00"

fly_start = "1979-06-01T00:00:00"
fly_end   = "1979-08-01T00:00:00"

arr_start = "1981-08-01T00:00:00"
arr_end   = "1981-10-01T00:00:00"

t0_d = PS.utc_to_et(dep_start); t1_d = PS.utc_to_et(dep_end)
t0_f = PS.utc_to_et(fly_start); t1_f = PS.utc_to_et(fly_end)
t0_a = PS.utc_to_et(arr_start); t1_a = PS.utc_to_et(arr_end)

dt = 1.0 * 86400.0 # Risoluzione 1 giorno

println("\n=== VOYAGER 2 RECREATION ===")

# ------------------------------------------------------------------
# 4. SOLVER
# ------------------------------------------------------------------
res = PS.porkchop_flyby(sys, bm_dep, bm_fly, bm_arr,
    (t0_d, t1_d, dt),
    (t0_f, t1_f, dt),
    (t0_a, t1_a, dt);
    R_fly_min = jupiter.radius + 50000.0,
    rpark_dep = 6678.0, 
    rpark_arr = NaN, 
    tof1_min  = 500 * 86400.0, tof1_max  = 800 * 86400.0,
    tof2_min  = 500 * 86400.0, tof2_max  = 1000 * 86400.0,
    dv_cap    = Inf
)

# ------------------------------------------------------------------
# 5. TROVA IL MIGLIORE (Launch + Flyby Cost only)
# ------------------------------------------------------------------
min_mission_dv = Inf
best_idx = (0,0)

for i in axes(res.dv_total, 1), k in axes(res.dv_total, 2)
    # Costo reale Voyager: Lancio + Correzione a Giove (Ignoriamo arrivo)
    cost = res.dv_dep[i,k] + res.dv_fly[i,k]
    
    if isfinite(cost) && cost < min_mission_dv
        min_mission_dv = cost
        best_idx = (i,k)
    end
end

if min_mission_dv == Inf
    println("[!] No valid trajectory found.")
    return
end

i, k = best_idx
t_dep_opt = res.tdep[i]
t_fly_opt = res.best_tfly[i,k]
t_arr_opt = res.tarr[k]

println("\n=== TRAIETTORIA OTTIMALE ===")
println("  Departure:   ", PS.et_to_utc(t_dep_opt; prec=0))
println("  Flyby:       ", PS.et_to_utc(t_fly_opt; prec=0))
println("  Arrival:     ", PS.et_to_utc(t_arr_opt; prec=0))
println("  Mission dV:  $(round(min_mission_dv, digits=4)) km/s")

# ------------------------------------------------------------------
# 6. PLOT 1: 2D PORKCHOP (Visualizzazione Costo Missione)
# ------------------------------------------------------------------
println("\nGenerazione Plot 2D...")
Z_mission = res.dv_dep .+ res.dv_fly # Matrice custom per il plot

fig_pork = PS.plot_grid(res;
    t0=t0_d, epoch0=Date(dep_start[1:10]),
    Z=Z_mission, 
    title="Voyager 2 (Launch + Flyby Cost)",
    zlabel="Effective dV (km/s)",
    zrange=(6.5, 9.0), 
    contour_levels=collect(6.5:0.25:9.0)
)

pork_path = joinpath(@__DIR__, "..", "plots", "voyager2_porkchop_2d.png")
save(pork_path, fig_pork)
println("Salvataggio 2D completato: $pork_path")

# ------------------------------------------------------------------
# 7. PLOT 2: 3D TRAJECTORY (Usando il TUO viz_traj.jl)
# ------------------------------------------------------------------
println("Generazione Plot 3D...")

# Chiamiamo la funzione che hai già in viz_traj.jl
# Assumiamo che PorkchopSolver esporti plot_trajectory_3d
fig_traj = PS.plot_trajectory_3d(
    sys, 
    (bm_dep, bm_fly, bm_arr), 
    (t_dep_opt, t_fly_opt, t_arr_opt);
    title_str = "Voyager 2: Earth -> Jupiter -> Saturn"
)

traj_path = joinpath(@__DIR__, "..", "plots", "voyager2_trajectory_3d.png")
save(traj_path, fig_traj)
println("Salvataggio 3D completato: $traj_path")
