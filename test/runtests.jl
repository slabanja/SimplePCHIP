using PCHIPInterpolation
using Test

using ForwardDiff: derivative
using Plots: plot

function test_interpolation_is_piecewise_monotone(xs, ys, N=10000)
    itp = @inferred Interpolator(xs, ys)
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


# Make sure the correct interval is identified
p = Interpolator([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
for search ∈ (PCHIPInterpolation._pchip_index_linear_search, PCHIPInterpolation._pchip_index_bisectional_search)
    @test @inferred search(p, 1.0) == 1
    @test @inferred search(p, 1.0 + eps(1.0)) == 1
    @test @inferred search(p, 2.0 - eps(2.0)) == 1
    @test @inferred search(p, 2.0) == 2
    @test @inferred search(p, 2.0 + eps(2.0)) == 2
    @test @inferred search(p, 3.0 - eps(3.0)) == 2
    @test @inferred search(p, 3.0) == 3
end


# Test interpolation
xs = [0.0,  1.2,  2.0,  5.0, 10.0, 11.0]
ys = [2.0,  2.1,  1.0,  0.0,  0.0,  3.0]
p = @inferred Interpolator(xs, ys)
@test @inferred p.(xs) == ys


# Test integration
p = @inferred Interpolator([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
@test @inferred integrate(p, 1, 3) == 0
@test @inferred integrate(p, 3, 1) == 0
@test @inferred integrate(p, 1, 2) == 0
@test @inferred integrate(p, 2, 3) == 0
@test @inferred integrate(p, 1.25, 2) == 0
@test @inferred integrate(p, 1, 2.25) == 0
@test @inferred integrate(p, 1.25, 2.25) == 0

p = @inferred Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
@test integrate(p, 1, 4) == 3*3/2 + 3
@test integrate(p, 1, 2) + integrate(p, 2, 4) == integrate(p, 1, 4)
@test integrate(p, 1, 2.75) + integrate(p, 2.75, 4) == integrate(p, 1, 4)
@test integrate(p, 3, 1) == -integrate(p, 1, 3)


# Test with ForwardDiff
p = @inferred Interpolator([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
@test @inferred derivative(p, 1) == 0
@test @inferred derivative(p, 2) == 0
@test @inferred derivative(p, 3) == 0
@test @inferred derivative(p, 1.25) == 0
@test @inferred derivative(p, 2.25) == 0

p = @inferred Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
@test @inferred derivative(p, 1) == -1
@test @inferred derivative(p, 2) == -1
@test @inferred derivative(p, 3) == -1
@test @inferred derivative(p, 4) == -1
@test @inferred derivative(p, 2.75) == -1

xs = [0.0,  1.2,  2.0,  5.0, 10.0, 11.0]
ys = [2.0,  2.1,  1.0,  0.0,  0.0,  3.0]
p = @inferred Interpolator(xs, ys)
for x ∈ [0, 1, 3.14, 5, 9, 11]
    @test derivative(b -> integrate(p, 0, b), x) ≈ p(x)
    @test derivative(b -> integrate(p, 1.1, b), x) ≈ p(x)
    @test derivative(b -> integrate(p, 11, b), x) ≈ p(x)
end


# Test out of domain
p = Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
@test p(1) == 4
@test p(4) == 1
@test_throws DomainError p(1 - 1e-6)
@test_throws DomainError p(4 + 1e-6)
@test_throws DomainError integrate(p, 1 - 1e-6, 4 + 1e-6)
@test_throws DomainError integrate(p, 1 - 1e-6, 3)
@test_throws DomainError integrate(p, 2, 4 + 1e-6)

# Test too short
@test_throws ArgumentError Interpolator([1.0], [4.0])
@test_throws ArgumentError Interpolator([], [])

# Test mismatched lengths
@test_throws DimensionMismatch Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0])
@test_throws DimensionMismatch Interpolator([1.0, 2.0, 3.0], [4.0, 3.0, 2.0, 1.0])

# Test non-increasing xs
@test_throws ArgumentError Interpolator([1.0, 2.0, 3.0, 3.0], [4.0, 3.0, 2.0, 1.0])


# Test plot recipe
xs = [0.0,  1.2,  2.0,  5.0, 10.0, 11.0]
ys = [2.0,  2.1,  1.0,  0.0,  0.0,  3.0]
p = @inferred Interpolator(xs, ys)
plot(p)
plot(p, markershape=:auto)


# Regression test for https://github.com/gerlero/PCHIPInterpolation.jl/issues/8
x = collect(range(0, stop=1, length=199))
itp = Interpolator(x, x)
@test itp(1) == 1
x = collect(range(0, stop=1, length=200))
itp = Interpolator(x, x)
@test itp(1) == 1
x = collect(range(0, stop=1, length=1000))
itp = Interpolator(x, x)
@test itp(1) == 1

# Regression test for https://github.com/gerlero/PCHIPInterpolation.jl/issues/11
v = [1, 2, 3]
p = Interpolator(v, v)
xs = copy(p.xs)
ys = copy(p.ys)
v .= -1
@test p.xs == xs
@test p.ys == ys
