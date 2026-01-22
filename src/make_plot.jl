# src/make_plots.jl
using CairoMakie
using JLD2
using Dates

# We must load the module so Julia understands what "PorkchopResult" is
include("PorkchopSolver.jl")
using .PorkchopSolver
PS = PorkchopSolver

println("Loading simulation data...")
data_path = joinpath(@__DIR__, "..", "mission_data.jld2")

if !isfile(data_path)
    error("Data file not found! Run demo.jl first.")
end

# Load variables back into memory
# quantities are loaded as a NamedTuple or Dictionary-like object
file = jldopen(data_path, "r")
res = file["res"]
tdep0 = file["tdep0"]
best_dv = file["best_dv"]
close(file)

println("Data loaded. Best dV found: $best_dv km/s")

# Setup Plotting Paths
plot_dir = joinpath(@__DIR__, "..", "plots")
mkpath(plot_dir)
epoch2026 = Date(2026, 1, 1)

# --- GENERATE PLOTS ---

println("1. Generating Porkchop (Total dV)...")
fig_dv = PS.plot_dv_total(res;
    t0 = tdep0, 
    epoch0 = epoch2026,
    levels = collect(5.0:0.5:15.0),
    zrange = (5.0, 15.0)
)
save(joinpath(plot_dir, "mars_2026_dv.png"), fig_dv)

println("2. Generating Launch C3...")
fig_c3 = PS.plot_c3(res;
    t0 = tdep0,
    epoch0 = epoch2026,
    levels = [10, 15, 20, 25, 30, 40, 50],
    zrange = (0.0, 50.0)
)
save(joinpath(plot_dir, "mars_2026_c3.png"), fig_c3)

println("3. Generating Arrival V-Inf...")
fig_vinf = PS.plot_vinf_arr(res;
    t0 = tdep0,
    epoch0 = epoch2026,
    levels = [2.0, 3.0, 4.0, 5.0, 6.0],
    zrange = (2.0, 8.0)
)
save(joinpath(plot_dir, "mars_2026_vinf.png"), fig_vinf)

println("Done! Plots saved to $plot_dir")