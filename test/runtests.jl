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
xs = 1.0:3.0
for xs in (xs, collect(xs)) 
    for search in (x -> (@inferred PCHIPInterpolation._findinterval_base(xs, x)),
                   x -> (@inferred PCHIPInterpolation._findinterval_custom(xs, x)))
        @test_throws DomainError search(0.0)
        @test_throws DomainError search(1.0 - eps(1.0))
        @test search(1.0) == 1
        @test search(1.0 + eps(1.0)) == 1
        @test search(2.0 - eps(2.0)) == 1
        @test search(2.0) == 2
        @test search(2.0 + eps(2.0)) == 2
        @test search(3.0 - eps(3.0)) == 2
        @test search(3.0) == 2
        @test_throws DomainError search(3.0 + eps(3.0))
        @test_throws DomainError search(4.0)
    end
end


# Test interpolation
xs = [0.0,  1.2,  2.0,  5.0, 10.0, 11.0]
ys = [2.0,  2.1,  1.0,  0.0,  0.0,  3.0]
itp = @inferred Interpolator(xs, ys)
@test itp.(xs) == ys
@inferred itp(xs[begin])


# Test interpolation with ranges
xs = 1:10
ys = 2:11
itp = @inferred Interpolator(xs, ys)
@test itp.(xs) == ys
@inferred itp(xs[begin])

xs = 2:101
ys = collect(3:102)
itp = @inferred Interpolator(xs, ys)
@test itp.(xs) == ys
@inferred itp(xs[begin])

itp2 = @inferred Interpolator(collect(xs), ys)

xs2 = [2.73, 3.14, 10, 15, 50, 60, 85]
@test itp2.(xs2) == itp.(xs2)
@inferred itp(xs2[begin])


# Test interpolator with custom derivatives (non-monotonic)
xs = 1:4
ys = xs
ds = -ones(length(ys))
itp = @inferred Interpolator(xs, ys, ds)
@test all(derivative(itp, x) == d for (x,d) in zip(xs,ds))
@test itp(1) < itp(1.9) > itp(2) < itp(2.9) > itp(3) < itp(3.9) > itp(4)


# Test integration
itp = @inferred Interpolator([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
@test (@inferred integrate(itp, 1, 3)) == 0
@test (@inferred integrate(itp, 3, 1)) == 0
@test (@inferred integrate(itp, 1, 2)) == 0
@test (@inferred integrate(itp, 2, 3)) == 0
@test (@inferred integrate(itp, 1.25, 2)) == 0
@test (@inferred integrate(itp, 1, 2.25)) == 0
@test (@inferred integrate(itp, 1.25, 2.25)) == 0

itp = @inferred Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
@test integrate(itp, 1, 4) == 3*3/2 + 3
@test integrate(itp, 1, 2) + integrate(itp, 2, 4) == integrate(itp, 1, 4)
@test integrate(itp, 1, 2.75) + integrate(itp, 2.75, 4) == integrate(itp, 1, 4)
@test integrate(itp, 3, 1) == -integrate(itp, 1, 3)


# Test with ForwardDiff
itp = @inferred Interpolator([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
@test (@inferred derivative(itp, 1)) == 0
@test (@inferred derivative(itp, 2)) == 0
@test (@inferred derivative(itp, 3)) == 0
@test (@inferred derivative(itp, 1.25)) == 0
@test (@inferred derivative(itp, 2.25)) == 0

itp = @inferred Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
@test (@inferred derivative(itp, 1)) == -1
@test (@inferred derivative(itp, 2)) == -1
@test (@inferred derivative(itp, 3)) == -1
@test (@inferred derivative(itp, 4)) == -1
@test (@inferred derivative(itp, 2.75)) == -1

xs = [0.0,  1.2,  2.0,  5.0, 10.0, 11.0]
ys = [2.0,  2.1,  1.0,  0.0,  0.0,  3.0]
itp = @inferred Interpolator(xs, ys)
for x ∈ [0, 1, 3.14, 5, 9, 11]
    @test derivative(b -> integrate(itp, 0, b), x) ≈ itp(x)
    @test derivative(b -> integrate(itp, 1.1, b), x) ≈ itp(x)
    @test derivative(b -> integrate(itp, 11, b), x) ≈ itp(x)
end


# Test out of domain
itp = Interpolator([1.0, 2.0, 3.0, 4.0], [4.0, 3.0, 2.0, 1.0])
@test itp(1) == 4
@test itp(4) == 1
@test_throws DomainError itp(1 - 1e-6)
@test_throws DomainError itp(4 + 1e-6)
@test_throws DomainError integrate(itp, 1 - 1e-6, 4 + 1e-6)
@test_throws DomainError integrate(itp, 1 - 1e-6, 3)
@test_throws DomainError integrate(itp, 2, 4 + 1e-6)

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
itp = @inferred Interpolator(xs, ys)
plot(itp)
plot(itp, markershape=:auto)


# Regression test for https://github.com/gerlero/PCHIPInterpolation.jl/issues/8
x = collect(range(0, stop=1, length=24))
itp = Interpolator(x, x)
@test itp(1) == 1
x = collect(range(0, stop=1, length=25))
itp = Interpolator(x, x)
@test itp(1) == 1
x = collect(range(0, stop=1, length=100))
itp = Interpolator(x, x)
@test itp(1) == 1

# Regression test for https://github.com/gerlero/PCHIPInterpolation.jl/issues/11
v = [1, 2, 3]
itp = Interpolator(v, v)
xs = copy(itp.xs)
ys = copy(itp.ys)
v .= -1
@test itp.xs == xs
@test itp.ys == ys

# Test correctness of _pchip_ds_scipy implementation following https://github.com/gerlero/PCHIPInterpolation.jl/issues/31
# Tests that the interpolated values are the same as SciPy
xs = [0.0,  1.2,  2.0,  5.0, 10.0, 11.0]
ys = [2.0,  2.1,  1.0,  0.0,  0.0,  3.0]

xps = [0.6, 1.6, 3.5, 7.5, 10.5]
yps = [2.08750000, 1.61081474, 0.27194471, 0.00000000, 1.06250000]

itp = @inferred Interpolator(xs, ys)
@test itp.(xps) ≈ yps atol=1e-8
