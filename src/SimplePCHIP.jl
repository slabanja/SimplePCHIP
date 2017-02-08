module SimplePCHIP
#
# Simple PCHIP implementation following Fritsch and Carlson, SIAM J. NUMER. ANAL. 17 (1980) 238-246.
# Derivative calculated in a fashion similar to SciPy's PchipInterpolate
#

immutable pchip
    N :: Integer
    xs :: Array{Float64}
    ys :: Array{Float64}
    ds :: Array{Float64}
end


function interp(pchip, xs)
    [_interp(pchip, x) for x ∈ xs]
end

ϕ(t) = 3t^2 - 2t^3
ψ(t) = t^3 - t^2

function _interp(pchip, x :: Number)
    i = _pchip_index(pchip, x)
    x1, x2 = pchip.xs[i:i+1]
    y1, y2 = pchip.ys[i:i+1]
    d1, d2 = pchip.ds[i:i+1]
    h = x2 - x1

    (y1 * ϕ((x2-x)/h)
     + y2 * ϕ((x-x1)/h)
     - d1*h * ψ((x2-x)/h)
     + d2*h * ψ((x-x1)/h))
end

function _pchip_index(pchip, x)
    N = pchip.N
    if N < 200  # Approximate performance cross-over on my old intel i7-3517U
        i = _pchip_index_linear_search(pchip, x)
    else
        i = _pchip_index_bisectional_search(pchip, x)
    end
    if i == N
        # Treat right endpoint as part of rightmost interval
        assert(x ≈ pchip.xs[N])
        i = N-1
    end
    i
end

function _pchip_index_linear_search(pchip, x)
    xmin = pchip.xs[1]
    assert(x >= xmin)

    i = 1
    N = pchip.N
    while i < N  &&  x >= pchip.xs[i+1]
        i = i + 1
    end
    i
end

function _pchip_index_bisectional_search(pchip, x)
    N = pchip.N
    imin, imax = 1, N
    xmin = pchip.xs[imin]
    xmax = pchip.xs[imax]
    assert(x >= xmin && x <= xmax)

    i = imin + div(imax - imin + 1, 2)
    while imin < imax
        if x < pchip.xs[i]
            imax = i - 1
        elseif x >= pchip.xs[i+1]
            imin = i + 1
        else
            break
        end
        i = imin + div(imax - imin + 1, 2)
    end
    i
end


"Similar to how SciPy's PCHIP does it"
function _initial_ds_scipy(xs, ys)
    h(i) = xs[i+1]-xs[i]
    Δ(i) = (ys[i+1]-ys[i]) / h(i)

    N = length(xs)
    ds = similar(xs)
    if N == 2
        ds[:] = Δ(1)
    else
        Δl = Δ(1)
        hl = h(1)
        for i ∈ 2:N-1
            Δr = Δ(i)
            hr = h(i)
            if sign(Δl) != sign(Δr) || Δl ≈ 0.0 || Δr ≈ 0.0  # isapprox(Δl, 0.0, atol=eps()) || isapprox(Δr, 0.0, atol=eps())
                ds[i] = 0.0
            else
                wl = 2hl + hr
                wr = hl + 2hr
                ds[i] = (wl + wr) / (wl/Δl + wr/Δr)
            end
            Δl = Δr
            hl = hr
        end
        ds[1] = _edge_derivative(h(1), h(2), Δ(1), Δ(2))
        ds[N] = _edge_derivative(h(N-1), h(N-2), Δ(N-1), Δ(N-2))
    end
    ds
end

function _edge_derivative(h1, h2, Δ1, Δ2)
    d = ((2h1 + h2)*Δ1 - h2*Δ2) / (h1 + h2)
    if sign(d) != sign(Δ1)
        d = 0.0
    elseif sign(Δ1) != sign(Δ2)  &&  abs(d) > abs(3Δ1)
        d = 3Δ1
    end
    d
end


function create_pchip(xs, ys)
    xs_ = [x for x ∈ xs]
    ys_ = [x for x ∈ ys]
    _assert_xs_ys(xs_, ys_)
    ds = _initial_ds_scipy(xs_, ys_)
    pchip(length(xs_), xs_, ys_, ds)
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
