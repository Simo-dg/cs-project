"""
    viz.jl

    This module provides 2D visualization utilities for porkchop plots and mission analysis in astrodynamics.
    It includes functions for plotting delta-v, C3, and v-infinity grids, as well as utilities for axis formatting and masking invalid data.
    The plots are generated using CairoMakie for high-quality output.

    Main features:
    - 2D grid plotting for porkchop analysis (plot_grid)
    - Axis and date formatting utilities
    - Data masking for invalid or infeasible solutions
"""


using CairoMakie
using Dates

@inline function _masked(Z::AbstractMatrix{<:Real}, ok::AbstractMatrix{Bool})
    M = Matrix{Float64}(undef, size(Z)...)
    @inbounds for i in axes(Z,1), j in axes(Z,2)
        v = Z[i,j]
        M[i,j] = (ok[i,j] && isfinite(v)) ? Float64(v) : NaN
    end
    return M
end

@inline function _axes_days(res, t0::Float64)
    dep_days = (res.tdep .- t0) ./ 86400.0
    arr_days = (res.tarr .- t0) ./ 86400.0
    return dep_days, arr_days
end

function _date_ticks(epoch0::Date, days::AbstractVector{<:Real}; nticks::Int=6)
    lo = minimum(days); hi = maximum(days)
    ticks = range(lo, hi; length=nticks)
    labels = string.(epoch0 .+ Day.(round.(Int, ticks)))
    return (collect(ticks), labels)
end

function plot_grid(res;
        t0::Float64,
        epoch0::Date,
        Z::AbstractMatrix,
        title::String,
        zlabel::String,
        zrange=nothing,
        contour_levels=nothing,
        contour_color=:black,
        contour_style=:solid,
        colormap=:viridis)  

    dep_days, arr_days = _axes_days(res, t0)
    ZZ = _masked(Z, res.ok)

    fig = Figure(size=(1100, 750))
    ax = Axis(fig[1,1],
        xlabel="Departure date",
        ylabel="Arrival date",
        title=title
    )

    hm = heatmap!(ax, dep_days, arr_days, ZZ;
        interpolate=false,
        nan_color=:transparent,
        colorrange = (zrange === nothing ? automatic : zrange),
        colormap = colormap  
    )
    Colorbar(fig[1,2], hm, label=zlabel)

    if contour_levels !== nothing
        contour!(ax, dep_days, arr_days, ZZ;
            levels=contour_levels,
            color=contour_color,
            linewidth=1.2,
            linestyle=contour_style
        )
    end

    ax.xticks = _date_ticks(epoch0, dep_days)
    ax.yticks = _date_ticks(epoch0, arr_days)

    fig
end

function plot_c3(res; t0, epoch0, levels=[10,20,30], zrange=nothing)
    plot_grid(res; t0=t0, epoch0=epoch0, Z=res.c3, title="Departure C3", zlabel="C3", zrange=zrange, contour_levels=levels)
end
function plot_vinf_arr(res; t0, epoch0, levels=[2,4,6], zrange=nothing)
    plot_grid(res; t0=t0, epoch0=epoch0, Z=res.vinf_arr, title="Arr v-inf", zlabel="km/s", zrange=zrange, contour_levels=levels)
end
function plot_dv_total(res; t0, epoch0, levels=[6,8,10], zrange=nothing)
    plot_grid(res; t0=t0, epoch0=epoch0, Z=res.dv, title="Total dV", zlabel="km/s", zrange=zrange, contour_levels=levels)
end