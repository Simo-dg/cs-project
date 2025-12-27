module PorkchopSolver

using LinearAlgebra
using StaticArrays
using Dates
using Printf

include("types.jl")
include("stumpff.jl")
include("kepler.jl")
include("lambert.jl")
include("ephemeris.jl")
include("porkchop.jl")
include("viz.jl")
include("gravity_assist.jl")

export Vec3, Body, BodyModel, KeplerianElements, TwoBodySystem, CircularCoplanarEphemeris, KeplerEphemeris
export PorkchopResult, porkchop_grid
export lambert_uv, propagate_universal
export plot_porkchop
export flyby_turning_angle, powered_flyby_delta_v
export porkchop_grid_dep_arr
export plot_c3, plot_vinf_arr, plot_dv_total
end # module