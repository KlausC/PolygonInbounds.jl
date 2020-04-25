# PolygonInbounds

[![Build Status](https://travis-ci.com/KlausC/PolygonInbounds.jl.svg?branch=master)](https://travis-ci.com/KlausC/PolygonInbounds.jl)
[![Codecov](https://codecov.io/gh/KlausC/PolygonInbounds.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/KlausC/PolygonInbounds.jl)

## INPOLY: Fast points-in-polygon query in Julia

Implementation and improvement of function INPOLY2 in Julia. Some new
features have been added.

The implementation claims to be fast for multiple points to check at once.
The cost is `O((N+M)*(log(M)+1))` where `M` is the number of points to check and `N`
the number of edges.

Link to original Matlab sources: [inpoly2.m](https://github.com/dengwirda/inpoly)
The original algorithm was developed by Darren Engwirda in 2017.
The Euclidian distance definition of "on-boundary" and the support for multiple
areas has been added.

## Description

```
    stat = inpoly2(points, nodes, edges[, atol=, rtol=])
```

determines for each point of `points`
two status bits: `inside` and `onboundary`. These bits are stored in the output
matrix `stat[:,1:2]` with the same row indices as `points`.

A point is considered "inside", if the ray starting from x heading east intersects an
odd number of edges.

A point is considered "on-boundary", if its Euclidian distance to any of the edges
of the polygon is less than `tol = max(atol, rtol*sizefactor)`. `sizefactor` is the
maximum extension of the smallest box containing all edges.

The cost factor depends severely on `tol`, if that is greater than the average
length of the edges.

`points` and `nodes` are matrices consisting of x, y in first and second column or
vectors of point-objects `p` with `p[1]` and `p[2]` the x and y coordinates.
`edges` is a matrix of indices into `nodes`, which define the egdges of a polygon
or collection of polygons.

The polygons may be unconnected and self-intersecting.

Each edge may be associated with one or more area indices, which are stored in adjacent
columns of `edges[:,3:end]`. If there is more than one additional area colum,
the output array becomes 3-dimensional with elements `stats[:,1:2,area]`. Here `area`
are the area indices as stored in the additional columns of `edges`.

The definitions of "inside" and "on-boundary" related to an area consider only edges,
which have this area associated in one of the added columns of `edges`.
Area index `0` indicates unused. It is possible to assign each edge to zero, one, or
more areas in this way.

The `polydemo` functions produce plots of some selected examples and may be used
to enjoy.

### Usage:

```julia

Pkg>add PolygonInbounds

using Plots
using PolygonInbounds
using .Demo
Demo.setplot(Plots)

polydemo(1; tol=0.2); # r is the number of points on a grid
polydemo(2; r = 10^5, tol=0.01); # polygons of the "lakes" mesh data set 
polydemo(3; r = 10^5, tol=0.01); # polygons of the "coast" mesh data set 
polydemo(4; tol=0.1) # simple tetragon 
polydemo(5; tol=0.2) # two polygons with different areas 

points = [0.05 0.0; 1 1; -1 1]
nodes =  [0.0 0; 0 10; 10 10; 10 0]
edges = [1 2; 2 3; 3 4; 4 1]
tol = 1e-1

stat = inpoly2(points, nodes, edges, atol=tol)

```

Docu of the original Matlab implementation: [inpoly2](./doc/original.md).

