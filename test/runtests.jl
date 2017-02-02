import SimplePCHIP
using Base.Test

reload("SimplePCHIP")
s = SimplePCHIP


function test_interpolation_is_piecewise_monotone(xs, ys, N=100)
    pchip = s.create_pchip(xs, ys)
    for (i,j) in monotone_intervals(ys)
        @test is_monotone(s.interp(pchip,
                                   linspace(xs[i], xs[j], N)))
    end
end

function monotone_intervals(xs)
    N = length(xs)
    function _monotone_intervals()
        up_down = 0
        previous_i = 1
        for i = 1:N-1
            current_x, next_x = xs[i:i+1]
            if up_down != 0  &&  up_down*current_x > up_down*next_x
                produce((previous_i, i))
                previous_i = i
                up_down = 0
            end
            if up_down == 0  &&  current_x != next_x
                up_down = current_x < next_x ? 1 : -1
            end
        end
        produce((previous_i, N))
    end

    [interval for interval in Task(_monotone_intervals)]
end

function is_monotone(zs)
    (all(is_approx_le(z, znext) for (z,znext) = zip(zs[1:end-1], zs[2:end]))
    || all(is_approx_le(znext, z) for (z,znext) = zip(zs[1:end-1], zs[2:end])))
end

function is_approx_le(a, b)
    a < b || isapprox(a, b, atol=eps())
end

xs = [0.0 1.0]
ys = [1.0 1.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [0.0 1.0]
ys = [1.0 0.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [0.0 1.0]
ys = [0.0 1.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [0.0 1.0 2.0 100.0]
ys = [0.0 0.0 1.0   1.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [1.0 2.0 5.0 6.0 10.0 22.0]
ys = [0.0 1.0 0.0 1.0  0.0  1.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [   0.0 10.0 10000.0]
ys = [-100.0  1.2     1.1]
test_interpolation_is_piecewise_monotone(xs, ys, N=10000)


# Example data from Fritsch and Carlson, SIAM J. NUMER. ANAL. 17 (1980) 238-246.
xs = [ 0.0  2  3  5  6  8  9   11 12 14 15]
ys = [10.0 10 10 10 10 10 10.5 15 50 60 85]
test_interpolation_is_piecewise_monotone(xs, ys)


# Test the test helper functions
@test is_monotone([0.0 0.0 0.0])
@test is_monotone([0.0 1.0 10.0])
@test is_monotone([10.0 1.0 0.0])
@test ! is_monotone([0.0 1.0 0.0])
@test ! is_monotone([1.0 0.0 1.0])
