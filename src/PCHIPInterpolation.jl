module PCHIPInterpolation

export Interpolator, integrate


_is_strictly_increasing(xs::AbstractVector) = all(a < b for (a,b) in zip(xs[begin:end-1], xs[begin+1:end]))
_is_strictly_increasing(xs::AbstractRange) = step(xs) > 0

function _pchip_edge_derivative(h1, h2, Δ1, Δ2)
    d = ((2h1 + h2)*Δ1 - h2*Δ2) / (h1 + h2)
    if sign(d) != sign(Δ1)
        d = 0.0
    elseif sign(Δ1) != sign(Δ2)  &&  abs(d) > abs(3Δ1)
        d = 3Δ1
    end
    d
end

function _pchip_ds_scipy(xs::AbstractVector, ys::AbstractVector)
    h(i) = xs[i+1]-xs[i]
    Δ(i) = (ys[i+1]-ys[i]) / h(i)

    N = length(xs)
    ds = similar(ys./xs)
    if N == 2
        ds[:] .= Δ(1)
    else
        Δl = Δ(1)
        hl = h(1)
        for i ∈ 2:N-1
            Δr = Δ(i)
            hr = h(i)
            if sign(Δl) != sign(Δr) || Δl ≈ 0.0 || Δr ≈ 0.0
                ds[i] = 0.0
            else
                wl = 2hr + hl
                wr = hr + 2hl
                ds[i] = (wl + wr) / (wl/Δl + wr/Δr)
            end
            Δl = Δr
            hl = hr
        end
        ds[1] = _pchip_edge_derivative(h(1), h(2), Δ(1), Δ(2))
        ds[N] = _pchip_edge_derivative(h(N-1), h(N-2), Δ(N-1), Δ(N-2))
    end
    ds
end

struct Interpolator{Xs,Ys,Ds}
    xs::Xs
    ys::Ys
    ds::Ds

    function Interpolator(xs::AbstractVector, ys::AbstractVector)
        length(xs) ≥ 2 || throw(ArgumentError("xs must have at least 2 elements"))
        _is_strictly_increasing(xs) || throw(ArgumentError("xs must be strictly increasing"))
        length(xs) == length(ys) || throw(DimensionMismatch("xs and ys must have the same length"))

        xs = deepcopy(xs)
        ys = deepcopy(ys)
        ds = _pchip_ds_scipy(xs, ys)

        new{typeof(xs),typeof(ys),typeof(ds)}(xs, ys, ds)
    end

    function Interpolator(xs::AbstractVector, ys::AbstractVector, _ds::AbstractVector)
        length(xs) ≥ 2 || throw(ArgumentError("xs must have at least 2 elements"))
        _is_strictly_increasing(xs) || throw(ArgumentError("xs must be strictly increasing"))
        length(xs) == length(ys) == length(_ds) || throw(DimensionMismatch("xs, ys and _ds must have the same length"))

        xs = deepcopy(xs)
        ys = deepcopy(ys)
        ds = deepcopy(_ds)

        new{typeof(xs),typeof(ys),typeof(ds)}(xs, ys, ds)
    end
end


@inline function _findinterval_base(xs, x) # Generic binary search from Julia Base
    @boundscheck length(xs) ≥ 2 || throw(BoundsError(xs))

    i = !isnan(x) ? searchsortedlast(xs, x) : lastindex(xs) # searchsortedlast fails with NaN in Julia 1.5

    if i < firstindex(xs)
        throw(DomainError(x, "Below interpolation range"))
    end

    if i == lastindex(xs)
        if x > @inbounds xs[i]
            throw(DomainError(x, "Above interpolation range"))
        end
        i -= 1 # Treat right endpoint as part of rightmost interval
    end

    return i
end

@inline function _findinterval_custom(xs, x) # SciPy-like binary search from SimplePCHIP
    @boundscheck length(xs) ≥ 2 || throw(BoundsError(xs))

    imin = firstindex(xs)

    if x < @inbounds xs[imin]
        throw(DomainError(x, "Below interpolation range"))
    end

    imax = lastindex(xs)
    xmax = @inbounds xs[imax]

    if x > xmax
        throw(DomainError(x, "Above interpolation range"))
    elseif x == xmax
        return imax - 1 # Treat right endpoint as part of rightmost interval
    end

    i = imin + (imax - imin + 1)÷2
    while imin < imax
        if x < @inbounds xs[i]
            imax = i - 1
        elseif x >= @inbounds xs[i+1]
            imin = i + 1
        else
            break
        end
        i = imin + (imax - imin + 1)÷2
    end

    return i
end

@inline _findinterval(itp::Interpolator{<:AbstractRange}, x) = _findinterval_base(itp.xs, x) # Base binary search has an overload for ranges
@inline _findinterval(itp::Interpolator, x) = _findinterval_custom(itp.xs, x) # Otherwise, the custom search is preferred


@inline _x(::Interpolator, x) = x
@inline _x(itp::Interpolator, x, _) = _x(itp, x)
Base.@propagate_inbounds _x(itp::Interpolator, ::Val{:begin}, i) = itp.xs[i]
Base.@propagate_inbounds _x(itp::Interpolator, ::Val{:end}, i) = itp.xs[i+1]

@inline _evaluate(itp::Interpolator, ::Val{:begin}, i) = itp.ys[i]
@inline _evaluate(itp::Interpolator, ::Val{:end}, i) = itp.ys[i+1]

@inline _derivative(itp::Interpolator, ::Val{:begin}, i) = itp.ds[i]
@inline _derivative(itp::Interpolator, ::Val{:end}, i) = itp.ds[i+1]

@inline _ϕ(t) = 3t^2 - 2t^3
@inline _ψ(t) = t^3 - t^2

Base.@propagate_inbounds function _evaluate(itp::Interpolator, x, i)
    x1 = _x(itp, Val(:begin), i)
    x2 = _x(itp, Val(:end), i)
    h = x2 - x1

    y1 = _evaluate(itp, Val(:begin), i)
    y2 = _evaluate(itp, Val(:end), i)

    d1 = _derivative(itp, Val(:begin), i)
    d2 = _derivative(itp, Val(:end), i)

    return (y1 * _ϕ((x2-x)/h)
           + y2 * _ϕ((x-x1)/h)
           - d1*h * _ψ((x2-x)/h)
           + d2*h * _ψ((x-x1)/h))
end

@inline _evaluate(itp::Interpolator, x) = @inbounds _evaluate(itp, x, _findinterval(itp, x))

@inline (itp::Interpolator)(x::Number) = _evaluate(itp, x)


Base.@propagate_inbounds function _integrate(itp::Interpolator, a, b, i)
    a_ = _x(itp, a, i)
    b_ = _x(itp, b, i)
    return (b_ - a_)/6*(_evaluate(itp, a, i) + 4*_evaluate(itp, (a_ + b_)/2, i) + _evaluate(itp, b, i)) # Simpson's rule
end

Base.@propagate_inbounds function _integrate(itp::Interpolator, a, b, i, j)
    if i == j
        return _integrate(itp, a, b, i)
    end

    integral = _integrate(itp, a, Val(:end), i)
    for k in i+1:j-1
        integral += _integrate(itp, Val(:begin), Val(:end), k)
    end
    integral += _integrate(itp, Val(:begin), b, j)

    return integral
end

@inline function _integrate(itp::Interpolator, a, b)
    if b < a
        return -_integrate(itp, b, a)
    end

    return @inbounds _integrate(itp, a, b, _findinterval(itp, a), _findinterval(itp, b))
end

@inline integrate(itp::Interpolator, a::Number, b::Number) = _integrate(itp, a, b)


if !isdefined(Base, :get_extension) # For compatibility with Julia < 1.9
    include("../ext/PCHIPInterpolationRecipesBaseExt.jl")
end

end
