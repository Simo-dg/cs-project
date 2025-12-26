using CairoMakie
using Dates

# Fill NaNs with +Inf for contour stability (they’ll just be “outside”)
@inline function _nan_to_inf!(A::AbstractMatrix{<:Real})
    @inbounds for i in eachindex(A)
        if !isfinite(A[i])
            A[i] = Inf
        end
    end
    return A
end

function plot_porkchop(res::PorkchopResult;
        epoch0::Date=Date(2026,1,1),
        title::String="Porkchop (Δv)",
        c3_levels = [0,5,10,15,20,30,40,60,80,100],
        vinf_levels = [1,2,3,4,5,6,7,8,9,10])

    dep_days = res.tdep ./ 86400.0
    arr_days = res.tarr ./ 86400.0

    # Makie expects z dims = (length(x), length(y))
    # Our matrices are (nt, na) => need transpose for plotting with (dep_days, arr_days)
    Zdv  = copy(res.dv')        # (na, nt)?? careful: dv is (nt, na), so dv' is (na, nt)
    Zc3  = copy(res.c3')
    Zvin = copy(res.vinf_arr')

    # After transpose, we want dims (length(dep), length(arr)) = (nt, na).
    # dv' is (na, nt), so we must NOT transpose that way.
    # Correct: keep original orientation and transpose at call site instead.
    # Let's rebuild correctly:

    Zdv  = copy(res.dv)         # (nt, na)
    Zc3  = copy(res.c3)
    Zvin = copy(res.vinf_arr)

    # Mask invalid points
    @inbounds for i in axes(Zdv,1), j in axes(Zdv,2)
        if !res.ok[i,j] || !isfinite(Zdv[i,j])
            Zdv[i,j]  = NaN
            Zc3[i,j]  = NaN
            Zvin[i,j] = NaN
        end
    end

    fig = Figure(size=(1150, 780))
    ax = Axis(fig[1,1],
        xlabel="Departure date",
        ylabel="Arrival date",
        title=title
    )

    # heatmap accepts NaNs (we make them transparent)
    hm = heatmap!(ax, dep_days, arr_days, Zdv;
    interpolate=false,
    nan_color=:transparent,
    colorrange=(5.0, 20.0)
)
    Colorbar(fig[1,2], hm, label="Δv (km/s)")

    # contours: Contour.jl hates NaNs -> replace them with Inf before contouring
    Zc3c  = _nan_to_inf!(copy(Zc3))
    ZvinC = _nan_to_inf!(copy(Zvin))

    contour!(ax, dep_days, arr_days, Zc3c; levels=c3_levels, linewidth=1.0)
    contour!(ax, dep_days, arr_days, ZvinC; levels=vinf_levels, linewidth=1.0, linestyle=:dash)

    dep_ticks = range(minimum(dep_days), maximum(dep_days), length=6)
    arr_ticks = range(minimum(arr_days), maximum(arr_days), length=6)

    ax.xticks = (collect(dep_ticks), string.(epoch0 .+ Day.(round.(Int, dep_ticks))))
    ax.yticks = (collect(arr_ticks), string.(epoch0 .+ Day.(round.(Int, arr_ticks))))

    return fig
end