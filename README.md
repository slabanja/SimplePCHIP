# PCHIPInterpolation.jl

[![Build Status](https://github.com/gerlero/PCHIPInterpolation.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/gerlero/PCHIPInterpolation.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/gerlero/PCHIPInterpolation.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/gerlero/PCHIPInterpolation.jl)

PCHIP (Piecewise Cubic Hermite Interpolating Polynomial) spline interpolation of arbitrarily spaced one-dimensional data in Julia. This package is a fork of [SimplePCHIP](https://github.com/slabanja/SimplePCHIP) with some extra features.

PCHIP interpolation preserves monotonicity (i.e., it will not over- or undershoot monotonic data points). See [this SciPy documentation page](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html) for more details.


## Summary

### Create a PCHIP interpolator

```jl
using PCHIPInterpolation # Will ask you to install the package if you don't have it

xs = [0.0,  1.2,  2.0,  5.0, 10.0, 11.0]
ys = [2.0,  2.1,  1.0,  0.0,  0.0,  3.0]

itp = Interpolator(xs, ys)
```

### Evaluate

```jl
y = itp(1.5) # At a single point
ys = itp.(xs) # At multiple points
```

Attempts to evaluate outside the interpolation range will throw a [`DomainError`](https://docs.julialang.org/en/v1/base/base/#Core.DomainError) (i.e., the interpolator will not perform extrapolation).

### Plot (with [Plots](https://github.com/JuliaPlots/Plots.jl))

```jl
using Plots

plot(itp, markers=true, label="PCHIP")
```

![Plot](example.png)

The monotonicity-preserving property of PCHIP interpolation can be clearly seen in the plot.

### Compute a definite integral

```jl
integral = integrate(itp, 1, 3) # Integral between 1 and 3
```

### Compute a derivative (with [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl))

```jl
using ForwardDiff

dydx = ForwardDiff.derivative(itp, 1.5)
```
