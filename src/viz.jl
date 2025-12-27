# src/viz.jl
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
        contour_style=:solid)

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
        colorrange = (zrange === nothing ? automatic : zrange)
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

function plot_c3(res;
        t0::Float64,
        epoch0::Date,
        levels=[5,10,15,20,25,30,40,60],
        zrange=nothing)

    plot_grid(res;
        t0=t0, epoch0=epoch0,
        Z=res.c3,
        title="Earth → Mars porkchop — Departure C3 (km²/s²)",
        zlabel="C3 (km²/s²)",
        zrange=zrange,
        contour_levels=levels,
        contour_color=:black,
        contour_style=:solid
    )
end

function plot_vinf_arr(res;
        t0::Float64,
        epoch0::Date,
        levels=[2,3,4,5,6,7,8,9,10],
        zrange=nothing)

    plot_grid(res;
        t0=t0, epoch0=epoch0,
        Z=res.vinf_arr,
        title="Earth → Mars porkchop — Arrival v∞ (km/s)",
        zlabel="v∞,arr (km/s)",
        zrange=zrange,
        contour_levels=levels,
        contour_color=:black,
        contour_style=:dash
    )
end

function plot_dv_total(res;
        t0::Float64,
        epoch0::Date,
        levels=[6,7,8,9,10,11,12,13,14,15],
        zrange=nothing)

    plot_grid(res;
        t0=t0, epoch0=epoch0,
        Z=res.dv,
        title="Earth → Mars porkchop — Total Δv (km/s) (patched conics)",
        zlabel="Δv total (km/s)",
        zrange=zrange,
        contour_levels=levels,
        contour_color=:black,
        contour_style=:solid
    )
end