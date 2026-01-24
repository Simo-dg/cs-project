"""
	PorkchopSolver.jl

	This module provides a comprehensive toolkit for interplanetary trajectory design, porkchop plot analysis, and gravity assist mission planning.
	It integrates ephemeris models, Lambert solvers, porkchop plot computation, flyby analysis, and advanced visualization tools for astrodynamics applications.

	Main features:
	- Data structures for celestial bodies, ephemerides, and mission systems
	- Lambert and universal variable solvers
	- Porkchop plot computation and visualization
	- Gravity assist (flyby) analysis and powered flyby calculations
	- 3D trajectory visualization utilities
"""


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
include("flyby.jl")
include("viz_traj.jl")

export Vec3, Body, BodyModel, KeplerianElements, TwoBodySystem, CircularCoplanarEphemeris, KeplerEphemeris
export PorkchopResult, porkchop_grid
export lambert_uv, propagate_universal
export plot_porkchop
export flyby_turning_angle, powered_flyby_delta_v
export porkchop_grid_dep_arr
export plot_c3, plot_vinf_arr, plot_dv_total
export porkchop_flyby, FlybyResult
export plot_trajectory_3d

end 