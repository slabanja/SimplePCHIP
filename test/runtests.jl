using PCHIPInterpolation
using Test


function test_interpolation_is_piecewise_monotone(xs, ys, N=10000)
    itp = Interpolator(xs, ys)
    for (i,j) in monotone_intervals(ys)
        @test is_monotone([itp(x) for x ∈ range(xs[i], stop=xs[j], length=N)])
    end
end

function monotone_intervals(xs)
    N = length(xs)
    function _monotone_intervals(ch :: Channel)
        up_down = 0
        previous_i = 1
        for i = 1:N-1
            current_x, next_x = xs[i:i+1]
            if up_down != 0  &&  up_down*current_x > up_down*next_x
                put!(ch, (previous_i, i))
                previous_i = i
                up_down = 0
            end
            if up_down == 0  &&  current_x != next_x
                up_down = current_x < next_x ? 1 : -1
            end
        end
        put!(ch, (previous_i, N))
    end

    [interval for interval in Channel(_monotone_intervals)]
end

function is_monotone(zs)
    (all(is_approx_le(z, znext) for (z,znext) = zip(zs[1:end-1], zs[2:end]))
    || all(is_approx_le(znext, z) for (z,znext) = zip(zs[1:end-1], zs[2:end])))
end

function is_approx_le(a, b)
    a < b || a ≈ b
end

xs = [0.0, 1.0]
ys = [1.0, 1.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [0.0, 1.0]
ys = [1.0, 0.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [0.0, 1.0]
ys = [0.0, 1.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [0.0, 1.0, 2.0, 100.0]
ys = [0.0, 0.0, 1.0,   1.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [1.0, 2.0, 5.0, 6.0, 10.0, 22.0]
ys = [0.0, 1.0, 0.0, 1.0,  0.0,  1.0]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = [   0.0, 10.0, 10000.0]
ys = [-100.0,  1.2,     1.1]
test_interpolation_is_piecewise_monotone(xs, ys)

xs = range(0, stop=2π, length=5000)
ys = sin.(xs)
test_interpolation_is_piecewise_monotone(collect(xs), ys)

# Example data from Fritsch and Carlson, SIAM J. NUMER. ANAL. 17 (1980) 238-246.
xs = [ 0.0,  2,  3,  5,  6,  8,  9,   11, 12, 14, 15]
ys = [10.0, 10, 10, 10, 10, 10, 10.5, 15, 50, 60, 85]
test_interpolation_is_piecewise_monotone(xs, ys)


# Test the test helper functions
@test is_monotone([0.0, 0.0, 0.0])
@test is_monotone([0.0, 1.0, 10.0])
@test is_monotone([10.0, 1.0, 0.0])
@test ! is_monotone([0.0, 1.0, 0.0])
@test ! is_monotone([1.0, 0.0, 1.0])


# Test internal functions

# Make sure the correct interval is identified
p = Interpolator([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
for search ∈ (PCHIPInterpolation._pchip_index_linear_search, PCHIPInterpolation._pchip_index_bisectional_search)
    @test search(p, 1.0) == 1
    @test search(p, 1.0 + eps(1.0)) == 1
    @test search(p, 2.0 - eps(2.0)) == 1
    @test search(p, 2.0) == 2
    @test search(p, 2.0 + eps(2.0)) == 2
    @test search(p, 3.0 - eps(3.0)) == 2
    @test search(p, 3.0) == 3
end

# Test integration
@test integrate(p, 1, 3) == 0
@test integrate(p, 3, 1) == 0
@test integrate(p, 1, 2) == 0
@test integrate(p, 2, 3) == 0
@test integrate(p, 1.25, 2) == 0
@test integrate(p, 1, 2.25) == 0
@test integrate(p, 1.25, 2.25) == 0

p = Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
@test integrate(p, 1, 4) == 3*3/2 + 3
@test integrate(p, 1, 2) + integrate(p, 2, 4) == integrate(p, 1, 4)
@test integrate(p, 1, 2.75) + integrate(p, 2.75, 4) == integrate(p, 1, 4)
@test integrate(p, 3, 1) == -integrate(p, 1, 3)

# Test out of domain
p = Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
@test p(1) == 4
@test p(4) == 1
@test_throws DomainError p(1 - 1e-6)
@test_throws DomainError p(4 + 1e-6)
@test_throws DomainError integrate(p, 1 - 1e-6, 4 + 1e-6)
@test_throws DomainError integrate(p, 1 - 1e-6, 3)
@test_throws DomainError integrate(p, 2, 4 + 1e-6)