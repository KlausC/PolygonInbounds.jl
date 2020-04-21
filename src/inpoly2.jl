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
function inpoly2(vert, node, edge=zeros(Int); atol::T=0.0, rtol::T=NaN, outformat=InOnBit) where T<:AbstractFloat

    rtol = !isnan(rtol) ? rtol : iszero(atol) ? eps(T)^0.85 : zero(T)
    poly = PolygonMesh(node, edge)
    points = PointsInbound(vert)
    npoints = length(points)
    nedges = edgecount(poly)

    vmin = minimum(points)
    vmax = maximum(points)
    pmin = minimum(poly)
    pmax = maximum(poly)

    lbar = sum(pmax - pmin)
    tol = max(abs(rtol * lbar), abs(atol))
    
    ac = areacount(poly)
    stat = ac > 1 ? falses(npoints,2,ac) : falses(npoints,2)
    # flip coordinates so expected efford is minimal

    dvert = vmax - vmin
    ix = dvert[1] < dvert[2] ? 1 : 2
    iyperm = sortperm(points, 3 - ix)

    inpoly2!(points, iyperm, poly, 1:nedges, ix, stat)
   
    sub = subproblems(points, iyperm, poly, ix, tol, [10; 10]) 
    n = 0
    for (iypermk, epermk) in sub
        n += 1
        println("$n: points: $(length(iypermk)) edges: $(length(epermk)) $epermk")
        inpoly2!(points, iyperm, poly, epermk, ix, tol, stat)
    end

    convertout(outformat, InOnBit, stat)
end

"""
    inpoly2_mat(vert, node, edge, fTolx, ftoly, stats)

INPOLY2_MAT the local m-code version of the crossing-number
test. Loop over edges; do a binary-search for the first ve-
rtex that intersects with the edge y-range; do crossing-nu-
mber comparisons; break when the local y-range is exceeded.
"""
function inpoly2!(points, iyperm, poly, eperm, ix::Integer, stat::S) where {N,T<:AbstractFloat,S<:AbstractArray{Bool,N}}

    nvrt = length(iyperm)   # number of points to be checked
    iy = 3 - ix
    
    #----------------------------------- loop over polygon edges
    for epos in eperm

        inod = edgeindex(poly, epos, 1)  # from
        jnod = edgeindex(poly, epos, 2)  # to
        # swap order of vertices
        if vertex(poly, inod, iy) > vertex(poly, jnod, iy)
            inod, jnod = jnod, inod
        end

        #------------------------------- calc. edge bounding-box
        xone = vertex(poly, inod, ix)
        ymin = yone = vertex(poly, inod, iy)
        xtwo = vertex(poly, jnod, ix)
        ymax = ytwo = vertex(poly, jnod, iy)

        xmin = min(xone, xtwo)
        xmax = max(xone, xtwo)

        ydel = ytwo - yone
        xdel = xtwo - xone

        # find top points[:,iy] < ymin by binary search
        ilow = searchfirst(points, iy, iyperm, ymin)
        #------------------------------- calc. edge-intersection
        # loop over all points with y ∈ [ymin,ymax)
        for jpos = ilow:nvrt
            jorig = iyperm[jpos]
            ypos = vertex(points, jorig, iy)
            ypos > ymax && break 
            xpos = vertex(points, jorig, ix)

            if xpos >= xmin
                if xpos <= xmax
                    #--------- inside extended bounding box of edge
                    mul1 = ydel * (xpos - xone)
                    mul2 = xdel * (ypos - yone)
                    if mul1 < mul2
                        #----- left of line && ypos exact to avoid multiple counting
                        flipio!(poly, stat, jorig, epos)
                    end
                end
            else # xpos < xmin - left of bounding box
                flipio!(poly, stat, jorig, epos)
            end
        end
    end
    stat
end
function inpoly2!(points, iyperm, poly, eperm, ix::Integer, veps::T, stat::S) where {N,T<:AbstractFloat,S<:AbstractArray{Bool,N}}

    nvrt = length(iyperm)   # number of points to be checked
    nedg = length(eperm)   # number of edges of the polygon mesh
    vepsx = vepsy = veps
    iy = 3 - ix

    #----------------------------------- loop over polygon edges
    for epos in eperm

        inod = edgeindex(poly, epos, 1)  # from
        jnod = edgeindex(poly, epos, 2)  # to
        # swap order of vertices
        if vertex(poly, inod, iy) > vertex(poly, jnod, iy)
            inod, jnod = jnod, inod
        end

        #------------------------------- calc. edge bounding-box
        xone = vertex(poly, inod, ix)
        yone = vertex(poly, inod, iy)
        xtwo = vertex(poly, jnod, ix)
        ytwo = vertex(poly, jnod, iy)

        xmin0 = min(xone, xtwo)
        xmax0 = max(xone, xtwo)
        xmin = xmin0 - vepsx
        xmax = xmax0 + vepsx
        ymin = yone - vepsy
        ymax = ytwo + vepsy

        ydel = ytwo - yone
        xdel = xtwo - xone
        xysq = xdel^2 + ydel^2
        feps = sqrt(xysq) * veps

        # find top points[:,iy] < ymin by binary search
        ilow = searchfirst(points, iy, iyperm, ymin)
        #------------------------------- calc. edge-intersection
        # loop over all points with y ∈ [ymin,ymax)
        for jpos = ilow:nvrt
            jorig = iyperm[jpos]
            ypos = vertex(points, jorig, iy)
            ypos > ymax && break 
            xpos = vertex(points, jorig, ix)

            if xmin <= xpos <= xmax
                #--------- inside extended bounding box of edge
                mul1 = ydel * (xpos - xone)
                mul2 = xdel * (ypos - yone)
                if abs(mul2 - mul1) <= feps
                    #------- distance from line through edge less veps
                    mul3 = xdel * (2xpos-xone-xtwo) + ydel * (2ypos-yone-ytwo)
                    if abs(mul3) <= xysq ||
                        hypot(xpos- xone, ypos - yone) <= veps ||
                        hypot(xpos- xtwo, ypos - ytwo) <= veps
                        # ---- round boundaries around endpoints of edge
                        setonbounds!(poly, stat, jorig, epos)
                    end
                end
            end
        end
    end
    stat
end

"""
    search lowest iy coordinate >= ymin
"""
function searchfirst(points::PointsInbound, iy::Integer, iyperm, ymin)
    ilow = 0
    iupp = length(points) + 1
    @inbounds while ilow < iupp - 1
        imid = ilow + (iupp-ilow) >>> 0x01
        if vertex(points, iyperm[imid], iy) < ymin
            ilow = imid
        else
            iupp = imid
        end
    end
    iupp
end

"""
    flipio!(poly, stat, p, epos)

Flip boolean value of in-out-status of point `p` for all areas associated with edge `epos`. 
"""
function flipio!(poly::PolygonMesh, stat::AbstractArray{Bool}, i::Integer, j::Integer)
    statop!((stat,j,a) -> begin stat[j,1,a] = !stat[j,1,a]; end, poly, stat, i, j)
end

"""
    setonbounds!(poly, stat, p, epos)

Set boolean value of on-bounds-status of point `p` for all areas associated with edge `epos`. 
"""
function setonbounds!(poly::PolygonMesh, stat::AbstractArray{Bool}, i::Integer, j::Integer)
    statop!((stat, j,a) -> begin stat[j,2,a] = true; end, poly, stat, i, j)
end
function statop!(f!::Function, poly::PolygonMesh{A}, stat::AbstractArray{Bool,N}, p::Integer, ed::Integer) where {A,N}
    for k in 1:max(A, 1)
        area = N == 3 ? areaindex(poly, ed, k) : 1
        if area > 0
            f!(stat, p, area)
        end
    end
end

struct Box{S}
    bmin::S
    bmax::S
    div
    Box(bmin::S, bmax::S, div) where S = new{S}(bmin, bmax, div)
end

"""
    subproblems(...)
Generate a set of smaller problems by tiling the domain into rectangular pieces.
Discard the rectangles, which neither contain points nor polygon edges.
"""
function subproblems(points::PointsInbound, iyperm, poly::PolygonMesh, ix::Integer, tol::AbstractFloat, div)

    iy = 3 - ix
    vmin = minimum(points)
    vmax = maximum(points)
    pmin = minimum(poly) .- tol
    pmax = maximum(poly) .+ tol
    bmin = max.(pmin, vmin)
    bmax = min.(pmax, vmax)
    box = Box(bmin[[ix,iy]], bmax[[ix,iy]], div[[ix,iy]])
    res = Dict{Tuple{Int,Int},NamedTuple{(:po,:ed),NTuple{2,Vector{Int}}}}()

    for epos in 1:edgecount(poly)
        inod = edgeindex(poly, epos, 1)  # from
        jnod = edgeindex(poly, epos, 2)  # to
        xone = vertex(poly, inod, ix)
        yone = vertex(poly, inod, iy)
        xtwo = vertex(poly, jnod, ix)
        ytwo = vertex(poly, jnod, iy)
        xmin, xmax = minmax(xone, xtwo)
        ymin, ymax = minmax(yone, ytwo)
        boxes = boxof(box, xmin - tol, ymin - tol, xmax + tol, ymax + tol)
        for k in boxes
            println("($xmin,$ymin),($xmax,$ymax) intersects $k")
            pushedge!(res, k, epos)
        end
    end

    for p in iyperm
        xpos = vertex(points, p, ix)
        ypos = vertex(points, p, iy)
        boxes = boxof(box, xpos, ypos)
        for k in boxes
            pushpoint!(res, k, p)
        end
    end
    values(res)
end

function pushedge!(res, k, epos)
    tup = get!(res, k) do
        (po = Int[], ed = Int[])
    end
    push!(tup.ed, epos)
    nothing
end

function pushpoint!(res, k, p)
    tup = get(res, k, nothing)
    if tup !== nothing
        push!(tup.po, p)
    end
    nothing
end

function boxof1(bmin, bmax, x, k)
    bmax <= bmin ? 1 : Int(floor((x - bmin) / (bmax - bmin) * k)) + 1
end
function clip(a, b, c)
    a < b ? 0 : a > c ? 0 : a
end

function boxof2(box::Box, x, y)
    kx = boxof1(box.bmin[1], box.bmax[1], x, box.div[1])
    ky = boxof1(box.bmin[2], box.bmax[2], y, box.div[2])
    kx, ky
end

function boxof(box::Box, x, y)
    kx, ky = boxof2(box, x, y)
    kx = clip(kx, 1, box.div[1])
    ky = clip(ky, 1, box.div[2])
    kx != 0 && ky != 0 ? [(kx, ky)] : Tuple{Int,Int}[]
end

function boxof(box::Box, x1, y1, x2, y2)
    kx1, ky1 = boxof2(box, x1, y1)
    kx2, ky2 = boxof2(box, x2, y2)
    kx1 = max(kx1, 1); kx2 = min(kx2, box.div[1])
    ky1 = max(ky1, 1); ky2 = min(ky2, box.div[2])
    ((kx, ky) for kx in kx1:kx2 for ky in ky1:ky2)
end


