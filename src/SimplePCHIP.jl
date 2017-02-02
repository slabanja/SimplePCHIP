module SimplePCHIP
#
# Simple PCHIP implementation following Fritsch and Carlson, SIAM J. NUMER. ANAL. 17 (1980) 238-246.
#
immutable pchip
    N :: Integer
    xs :: Array{Float64}
    ys :: Array{Float64}
    ds :: Array{Float64}
end


function interp(pchip, xs)
    [_interp(pchip, x) for x = xs]
end

function _interp(pchip, x :: Number)
    i = _pchip_interval_index(pchip, x)
    xi, xip1 = pchip.xs[i:i+1]
    yi, yip1 = pchip.ys[i:i+1]
    di, dip1 = pchip.ds[i:i+1]
    hi = xip1 - xi

    (yi*ϕ((xip1-x)/hi)
     + yip1*ϕ((x-xi)/hi)
     - di*hi*ψ((xip1-x)/hi)
     + dip1*hi*ψ((x-xi)/hi))
end

ϕ(t) = 3t^2 - 2t^3
ψ(t) = t^3 - t^2

function _pchip_interval_index_linear_search(pchip, x)
    xmin = pchip.xs[1]
    assert(x >= xmin)

    i = 1
    N = pchip.N
    while i < N  &&  x > pchip.xs[i+1]
        i = i + 1
    end
    assert(i < N || x == pchip.xs[N])
    i
end

function _pchip_interval_index_bisectional_search(pchip, x)
    N = pchip.N
    xmin = pchip.xs[1]
    xmax = pchip.xs[N]
    assert(x >= xmin && x <= xmax)

    N = pchip.N
    i = div(N+1, 2)

    while x < pchip.xs[i]  ||  x >= pchip.xs[i+1]
        if x < pchip.xs[i]
            i = div(i+1, 2)
        elseif x > pchip.xs[i+1]
            i = i + div(N-i+1, 2)
        end
    end
    i
end

_pchip_interval_index = _pchip_interval_index_linear_search


function _create_initial_ds(xs, ys)
    Δ(i) = (ys[i+1]-ys[i]) / (xs[i+1]-xs[i])
    d3p(i, x) = _three_point_interpolate_derivative(x, xs[i-1], ys[i-1], xs[i], ys[i], xs[i+1], ys[i+1])

    N = length(xs)
    ds = similar(xs)
    if N == 2
        ds[:] = Δ(1)
    else
        ds[1] = d3p(2, xs[1])
        for i = 2:N-1
            ds[i] = d3p(i, xs[i])
        end
        ds[N] = d3p(N-1, xs[N])

        for i = 1:N-1
            delta = Δ(i)
            if isapprox(delta, 0.0, atol=eps())
                ds[i] = ds[i+1] = 0.0
            else
                if sign(ds[i]) == -sign(delta)
                    ds[i] = 0.0
                end
                if sign(ds[i+1]) == -sign(delta)
                    ds[i+1] = 0.0
                end
            end
        end
    end
    ds
end
function _three_point_interpolate_derivative(x, x1, y1, x2, y2, x3, y3)
    (y1*(2x-x2-x3)/((x1-x2)*(x1-x3))
    + y2*(2x-x3-x1)/((x2-x3)*(x2-x1))
    + y3*(2x-x1-x2)/((x3-x1)*(x3-x2)))
end

function _map_ds!(xs, ys, ds, d_map)
    Δ(i) = (ys[i+1]-ys[i]) / (xs[i+1]-xs[i])
    N = length(xs)
    for i = 1:N-1
        if ! isapprox(Δ(i), 0.0, atol=eps())
            α = ds[i] / Δ(i)
            β = ds[i+1] / Δ(i)
            τ = d_map(α, β)
            ds[i] = τ * ds[i]
            ds[i+1] = τ * ds[i+1]
        end
    end
end

function _d_map_0(α, β)
    0.0
end
function _d_map_s1(α, β)
    if α > 3 || β > 3
        3.0 / max(α, β)
    else
        1.0
    end
end
function _d_map_s2(α, β)
    if α^2 + β^2 > 9
        3/sqrt(α^2 + β^2)
    else
        1.0
    end
end
function _d_map_s3(α, β)
    if α + β > 3
        3.0/(α + β)
    else
        1.0
    end
end
function _d_map_s4(α, β)
    if 2α + β > 3 && α + 2β > 3
        if α > β
            3.0/(α + 2β)
        else
            3.0/(2α + β)
        end
    else
        1.0
    end
end

_derivative_mappers = Dict(
    "0" => _d_map_0,
    "s1" => _d_map_s1,
    "s2" => _d_map_s2,
    "s3" => _d_map_s3,
    "s4" => _d_map_s4,
    )

function create_pchip(_xs, _ys, _derivative_mapper="s3")
    xs = [x for x=_xs]
    ys = [x for x=_ys]
    _assert_xs_ys(xs, ys)
    ds = _create_initial_ds(xs, ys)
    _map_ds!(xs, ys, ds, _derivative_mappers[_derivative_mapper])
    pchip(length(xs), xs, ys, ds)
end

function _assert_xs_ys(xs, ys)
    N = length(xs)
    assert(N > 1)
    assert(N == length(ys))
    assert_monotonic_increase(xs)
end

function assert_monotonic_increase(xs)
    foldl((a,b) -> (assert(a < b); b), xs)
end


end  # module
