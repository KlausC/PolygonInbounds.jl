using Test
using PolygonInbounds

@test polydemo(4, r=10^5, tol=1e-2) === nothing
