# PCHIPInterpolation.jl

[![Build Status](https://github.com/gerlero/PCHIPInterpolation.jl/workflows/CI/badge.svg)](https://github.com/gerlero/PCHIPInterpolation.jl/actions)
[![Coverage](https://codecov.io/gh/gerlero/PCHIPInterpolation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gerlero/PCHIPInterpolation.jl)

This package provides functionality to perform piecewise cubic hermite interpolating polynomial (PCHIP)
interpolation of arbitrarily spaced 1-dimensional data.


## Summary
The basic use of PCHIPInterpolation can be illustrated with the following snippet:
```jl
using PCHIPInterpolation
using Gadfly

xs = [0.0  1.2  2.0  5.0 10.0 11.0]
ys = [2.0  2.1  1.0  0.0  0.0  3.0]
itp = Interpolator(xs, ys)

xrange = range(xs[1], stop=xs[end], length=100)
yinterpolated = [itp(x) for x ∈ xrange]

plot(layer(x=xrange, y=yinterpolated, Geom.line),
     layer(x=xs, y=ys, Geom.point))
```

![](https://cloud.githubusercontent.com/assets/154866/23104112/e94d7eda-f6c7-11e6-9108-888555ed8d6a.png)

Besides being evaluated, `Interpolator` objects can also be integrated. For instance,
```jl
integral = integrate(itp, 1, 3)
```
computes the definite integral of `itp` between 1 and 3.

## Why PCHIP?
PCHIP interpolation preserves monotonicity.
E.g., if input data points are monotonically increasing, so will the interpolated points.
Also, the interpolated points will not overshoot.

It can be illustrated by zooming in on the plot above, at around x=1.2 and x=10.0,

![x=1.2](https://cloud.githubusercontent.com/assets/154866/23104705/51a5ea66-f6d3-11e6-816f-4f16057428d3.png)
![x=10.0](https://cloud.githubusercontent.com/assets/154866/23104707/577107fa-f6d3-11e6-8832-c25b9a033ba3.png)


## See also

PCHIPInterpolation is a fork of [SimplePCHIP](https://https://github.com/slabanja/SimplePCHIP), a package created to provide interpolation similar to SciPy's
[PchipInterpolation](http://scipy.github.io/devdocs/generated/scipy.interpolate.PchipInterpolator.html).

For further details on PCHIP interpolation, there is of course a wikipedia article about [Monotone cubic interpolation](https://en.wikipedia.org/wiki/Monotone_cubic_interpolation), also [this pdf about interpolation](https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/moler/interp.pdf) (with a focus on Matlab) provides details on PCHIP interpolation.
