@inline function stumpff_C(z::Float64)
    az = abs(z)
    if az < 1e-8
        z2 = z*z
        return 0.5 - z/24 + z2/720 - z2*z/40320
    elseif z > 0
        s = sqrt(z)
        return (1 - cos(s))/z
    else
        s = sqrt(-z)
        return (cosh(s) - 1)/(-z)
    end
end

@inline function stumpff_S(z::Float64)
    az = abs(z)
    if az < 1e-8
        z2 = z*z
        return 1/6 - z/120 + z2/5040 - z2*z/362880
    elseif z > 0
        s = sqrt(z)
        return (s - sin(s))/(s^3)
    else
        s = sqrt(-z)
        return (sinh(s) - s)/(s^3)
    end
end