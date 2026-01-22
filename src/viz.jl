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
        contour_style=:solid,
        colormap=:viridis)  # <--- NUOVO ARGOMENTO AGGIUNTO QUI

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
        colormap = colormap  # <--- PASSAGGIO ALLA HEATMAP
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

# Wrapper per compatibilità (non strettamente necessari se si usa plot_grid direttamente)
function plot_c3(res; t0, epoch0, levels=[10,20,30], zrange=nothing)
    plot_grid(res; t0=t0, epoch0=epoch0, Z=res.c3, title="Departure C3", zlabel="C3", zrange=zrange, contour_levels=levels)
end
function plot_vinf_arr(res; t0, epoch0, levels=[2,4,6], zrange=nothing)
    plot_grid(res; t0=t0, epoch0=epoch0, Z=res.vinf_arr, title="Arr v-inf", zlabel="km/s", zrange=zrange, contour_levels=levels)
end
function plot_dv_total(res; t0, epoch0, levels=[6,8,10], zrange=nothing)
    plot_grid(res; t0=t0, epoch0=epoch0, Z=res.dv, title="Total dV", zlabel="km/s", zrange=zrange, contour_levels=levels)
end