# PolygonInbounds

[![Build Status](https://travis-ci.com/KlausC/PolygonInbounds.jl.svg?branch=master)](https://travis-ci.com/KlausC/PolygonInbounds.jl)
[![Codecov](https://codecov.io/gh/KlausC/PolygonInbounds.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/KlausC/PolygonInbounds.jl)

### Purpose

Only purpose of the is to implement and improve function INPOLY2 in Julia.

The implementation claims to be fast for multiple points to check at once.
The algorithm was developed by Darren Engwirda in 2017.

`stats = inpoly2(points, nodes, edges[, atol=, rtol=])` determines for each point in `points`
a status of

* -1 outside polygon
* +1 inside polygon 
* 0  on boundary of polygon

A point is considered on boundary, if its Euclidian distance to any of the edges
of the polygon is less than `tol = max(atol, rtol*sizefactor)`.

`points` and `nodes` are matrices consisting of x, y in first and second column.
`edges` is a matrix of indices into `nodes`, which define the egdges of the polygon.
The polygon may be unconnected and self-overlapping.

Link to original Matlab sources: [inpoly2.m](https://github.com/dengwirda/inpoly)


### Usage:

```julia

]add //https://github.com/KlausC/PolygonInbounds.jl

using PolygonInbounds

polydemo.(1:4)

points = [0.05 0.0; 1 1; -1 1]
nodes =  [0.0 0; 0 10; 10 10; 10 0]
edges = [1 2; 2 3; 3 4; 4 1]
tol = 1e-1

stat = inpoly2(points, nodes, edges, atol=tol)

```

It follows the docu from the original Matlab implementation.

```
INPOLY2 compute "points-in-polygon" queries.
    [STAT] = INPOLY2(VERT,NODE,EDGE) returns the "inside/ou-
    tside" status for a set of vertices VERT and a polygon
    {NODE,EDGE} embedded in a two-dimensional plane. General
    non-convex and multiply-connected polygonal regions can
    be handled. VERT is an N-by-2 array of XY coordinates to
    be tested. STAT is an associated N-by-1 logical array,
    with STAT(II) = TRUE if VERT(II,:) is an interior point.
    The polygonal region is defined as a piecewise-straight-
    line-graph, where NODE is an M-by-2 array of polygon ve-
    rtices and EDGE is a P-by-2 array of edge indexing. Each
    row in EDGE represents an edge of the polygon, such that
    NODE(EDGE(KK,1),:) and NODE(EDGE(KK,2),:) are the coord-
    inates of the endpoints of the KK-TH edge. If the argum-
    ent EDGE is omitted it assumed that the vertices in NODE
    are connected in ascending order.
 
    [STAT,BNDS] = INPOLY2(..., FTOL) also returns an N-by-1
    logical array BNDS, with BNDS(II) = TRUE if VERT(II,:)
    lies "on" a boundary segment, where FTOL is a floating-
    point tolerance for boundary comparisons. By default,
    FTOL = EPS ^ 0.85.
 
    See also INPOLYGON

    This algorithm is based on a "crossing-number" test, co-
    unting the number of times a line extending from each
    point past the right-most region of the polygon interse-
    cts with the polygonal boundary. Points with odd counts
    are "inside". A simple implementation requires that each
    edge intersection be checked for each point, leading to
    O(N*M) complexity...
 
    This implementation seeks to improve these bounds:
 
  * Sorting the query points by y-value and determining can-
    didate edge intersection sets via binary-search. Given a
    configuration with N test points, M edges and an average
    point-edge "overlap" of H, the overall complexity scales
    like O(M*H + M*LOG(N) + N*LOG(N)), where O(N*LOG(N))
    operations are required for sorting, O(M*LOG(N)) operat-
    ions required for the set of binary-searches, and O(M*H)
    operations required for the intersection tests, where H
    is typically small on average, such that H << N.
 
  * Carefully checking points against the bounding-box asso-
    ciated with each polygon edge. This minimises the number
    of calls to the (relatively) expensive edge intersection
    test.

    Darren Engwirda : 2017 --
    Email           : de2363@columbia.edu
    Last updated    : 27/10/2018
```