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

function _pchip_index_linear_search(pchip, x)
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

function _pchip_index_bisectional_search(pchip, x)
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

_pchip_index = _pchip_index_linear_search

#
# function _initial_ds(xs, ys)
#     Δ(i) = (ys[i+1]-ys[i]) / (xs[i+1]-xs[i])
#     d3p(i, x) = _three_point_interpolate_derivative(x, xs[i-1], ys[i-1], xs[i], ys[i], xs[i+1], ys[i+1])
#
#     N = length(xs)
#     ds = similar(xs)
#     if N == 2
#         ds[:] = Δ(1)
#     else
#         ds[1] = d3p(2, xs[1])
#         for i = 2:N-1
#             ds[i] = d3p(i, xs[i])
#         end
#         ds[N] = d3p(N-1, xs[N])
#
#         for i ∈ 1:N-1
#             delta = Δ(i)
#             if isapprox(delta, 0.0, atol=eps())
#                 ds[i] = ds[i+1] = 0.0
#             else
#                 if sign(ds[i]) == -sign(delta)
#                     ds[i] = 0.0
#                 end
#                 if sign(ds[i+1]) == -sign(delta)
#                     ds[i+1] = 0.0
#                 end
#             end
#         end
#     end
#     ds
# end
# function _three_point_interpolate_derivative(x, x1, y1, x2, y2, x3, y3)
#     (y1*(2x-x2-x3)/((x1-x2)*(x1-x3))
#     + y2*(2x-x3-x1)/((x2-x3)*(x2-x1))
#     + y3*(2x-x1-x2)/((x3-x1)*(x3-x2)))
# end

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
            if sign(Δl) != sign(Δr) || isapprox(Δl, 0.0, atol=eps()) || isapprox(Δr, 0.0, atol=eps())
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

# function _map_ds!(xs, ys, ds, d_map)
#     Δ(i) = (ys[i+1]-ys[i]) / (xs[i+1]-xs[i])
#     N = length(xs)
#     for i = 1:N-1
#         if ! isapprox(Δ(i), 0.0, atol=eps())
#             α = ds[i] / Δ(i)
#             β = ds[i+1] / Δ(i)
#             τ = d_map(α, β)
#             ds[i] = τ * ds[i]
#             ds[i+1] = τ * ds[i+1]
#         end
#     end
# end
#
# _d_scale_zero(α, β) = 0.0
# _d_scale_unit(α, β) = 1.0
# function _d_scale_s1(α, β)
#     if α > 3 || β > 3
#         3.0 / max(α, β)
#     else
#         1.0
#     end
# end
# function _d_scale_s2(α, β)
#     if α^2 + β^2 > 9
#         3 / sqrt(α^2 + β^2)
#     else
#         1.0
#     end
# end
# function _d_scale_s3(α, β)
#     if α + β > 3
#         3.0 / (α + β)
#     else
#         1.0
#     end
# end
# function _d_scale_s4(α, β)
#     if 2α + β > 3  &&  α + 2β > 3
#         if α > β
#             3.0 / (α + 2β)
#         else
#             3.0 / (2α + β)
#         end
#     else
#         1.0
#     end
# end
#
# _derivative_scalers = Dict(
#     "0" => _d_scale_zero,
#     "1" => _d_scale_unit,
#     "s1" => _d_scale_s1,
#     "s2" => _d_scale_s2,
#     "s3" => _d_scale_s3,
#     "s4" => _d_scale_s4,
#     )

function create_pchip(xs, ys)
    xs_ = [x for x ∈ xs]
    ys_ = [x for x ∈ ys]
    _assert_xs_ys(xs_, ys_)
    ds = _initial_ds_scipy(xs_, ys_)
    # ds = _initial_ds(xs_, ys_)
    # _scale_ds!(xs_, ys_, ds, _derivative_scalers[derivative_scaler])
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
