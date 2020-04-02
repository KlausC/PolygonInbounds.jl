"""
    s = inpoly2(vert, node, [edge=Nothing, atol=0.0, rtol=eps(), outformat=InOnOut{1,0,-1}])

Check all points defined by `vert` are inside, outside or on bounds of
polygon defined by `node` and `edge`.

`vert` and `node` are matrices with the x-coordinates in the first
and the y-coordinates in the second column.

`edge` is an integer matrix of same size as `node`. The first column contains
the index of the starting node, the second column the index of the ending node.
The polygon needs to be closed. It may contain unconnected cycles.

A point is considered on boundary, if its Euclidian distance to any edge is less
than or equal `max(atol, rtol * span)`, where span is the maximum extension of the
polygon in either direction.

If a point is considered on bounds, no further check is made to determine it it is inside or 
outside numerically.

Expected effort of the algorithm is `(N+M)*log(M)`, where `M` is the number of points
and `N` is the number of polygon edges.
"""
function inpoly2(vert, node, edge=zeros(Int); atol::T=0.0, rtol::T=NaN, outformat=InOnOut{1,0,-1}) where T<:AbstractFloat

    rtol = isnan(rtol) ? eps(T)^0.85 : rtol
    polygon = PolygonType(node, edge)
    points = PointsInbound(vert)
    nnod = length(node)
    nvrt = length(points)

    vmin = minimum(points)
    vmax = maximum(points)
    ddxy = vmax - vmin
    lbar = sum(ddxy) / 2
    # flip coordinates so y-span of points >= x-span of points
    flip = ddxy[1] > ddxy[2]
    ivec = sortperm(points, 2-flip)
    stat = view(fill(Int8(-1), nvrt), ivec) # -1 outside, 0 onbound, +1 inside
    statv = view(stat, ivec)

    tol = max(abs(rtol * lbar), abs(atol))
    inpoly2!(points, ivec, polygon, flip, tol, statv)
    convertout(outformat, InOnOut{1,0,-1}, stat)
end

"""
    inpoly2_mat(vert, node, edge, fTol, stats)

INPOLY2_MAT the local m-code version of the crossing-number
test. Loop over edges; do a binary-search for the first ve-
rtex that intersects with the edge y-range; do crossing-nu-
mber comparisons; break when the local y-range is exceeded.
"""
function inpoly2!(points, ivec, polygon, flip::Bool, veps::AbstractFloat, stat::AbstractVector{<:Integer})
    nvrt = length(points) # the points to be checked
    nnod = length(polygon) # the vertices of the polygon
    
    ix = flip + 1
    iy = 2 - flip
    #----------------------------------- loop over polygon edges
    for epos = 1:nnod

        inod = edge(polygon, epos, 1)  # from
        jnod = edge(polygon, epos, 2)  # to
        if vertex(polygon, inod, iy) > vertex(polygon, jnod, iy)
            inod, jnod = jnod, inod
        end

        #------------------------------- calc. edge bounding-box
        xone = vertex(polygon, inod, ix)
        yone = vertex(polygon, inod, iy)
        xtwo = vertex(polygon, jnod, ix)
        ytwo = vertex(polygon, jnod, iy)

        xmin = min(xone, xtwo) - veps
        xmax = max(xone, xtwo) + veps
        ymin = yone - veps
        ymax = ytwo + veps

        ydel = ytwo - yone
        xdel = xtwo - xone
        feps = veps * hypot(xdel, ydel)

        # find top points[:,iy] < ymin by binary search
        ilow = 1
        iupp = nvrt
        while ilow < iupp - 1
            imid = ilow + (iupp-ilow) ÷ 2
            if vertex(points, ivec[imid], iy) < ymin
                ilow = imid
            else
                iupp = imid
            end
        end
        while ilow > 0 && vertex(points, ivec[ilow], iy) >= ymin
            ilow = ilow - 1
        end

        #------------------------------- calc. edge-intersection
        # loop over all points with y ∈ [ymin,ymax]
        for jpos = ilow+1:nvrt
            stat[jpos] == 0 && continue
            ypos = vertex(points, ivec[jpos], iy)
            ypos > ymax && break 
            xpos = vertex(points, ivec[jpos], ix)

            if xpos >= xmin
                if xpos <= xmax
                    #--------- inside extended bounding box of edge
                    mul1 = ydel * (xpos - xone)
                    mul2 = xdel * (ypos - yone)
                    if abs(mul2 - mul1) <= feps
                        #------- distance from line through edge less veps
                        if !(xdel * (xpos - xone) < ydel * (yone - ypos) &&
                             hypot(xpos- xone, ypos - yone) > veps ||
                             xdel * (xpos - xtwo) > ydel * (ytwo - ypos) &&
                             hypot(xpos- xtwo, ypos - ytwo) > veps)
                            # ---- round boundaries around endpoints of edge
                            stat[jpos]= 0
                        elseif mul1 < mul2 && yone <= ypos < ytwo
                            #----- left of line && ypos exact to avoid multiple counting
                            stat[jpos] = -stat[jpos]
                        end
                    elseif mul1 < mul2 && yone <= ypos < ytwo
                        #----- left of line && ypos exact to avoid multiple counting
                        stat[jpos] = -stat[jpos]
                    end
                end
            else # xpos < xmin - left of bounding box
                if yone <= ypos <  ytwo
                    #----- ypos exact to avoid multiple counting
                    stat[jpos] = -stat[jpos]
                end
            end
        end
    end
    stat
end
