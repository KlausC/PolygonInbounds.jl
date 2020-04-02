"""
    function inpoly2(vert, node, edge[, tol])

Check all points defined by `vert` are inside, outside or on bounds of
polygon defined by `node` and `edge`.
`vert` and `node` are matrices with the x-coordinates in the first
and the y-coordinates in the second column.
`edge` is an integer matrix of same size a `node`. The first column contains
the index of the starting node, the second column the index of the ending node.
The polygon needs to be closed. It may contain unconnected cycles.
"""
function inpoly2(vert::AbstractMatrix{T},
                 node::AbstractMatrix{T},
                 edge::AbstractArray{<:Integer,N}=zeros(Int);
                 atol::T=0.0, rtol::T=eps(T)^0.85) where {T<:Real,N}

    nnod = size(node,1)
    nvrt = size(vert,1)
    if N == 0
        edge = [1:nnod-1 2:nnod; nnod 1]
    elseif N == 1
        edge = [1:length(edge) edge]
    end
    checkfor2(v, s) = size(v,2) == 2 || throw(ArgumentError("array $s needs 2 columns"))
    checkfor2(node, "node"); checkfor2(edge, "edge"); checkfor2(vert, "vert")
    if size(edge, 1) != nnod
        throw(ArgumentError("edges need same length as nodes"))
    end
    if N == 2 && !isperm(view(edge,:,1)) || N >= 1 && !isperm(view(edge,:,2))
        throw(ArgumentError("edges are not permutation"))
    end

    vmin = minimum(vert, dims=1)
    vmax = maximum(vert, dims=1)
    ddxy = vmax - vmin
    lbar = sum(ddxy) / 2
    # flip coordinates so y-span of points >= x-span of points
    flip = ddxy[1] > ddxy[2]
    ivec = sortperm(view(vert,:,2-flip))
    vert = view(vert,ivec,:)
    stat = view(fill(Int8(-1), nvrt), ivec) # -1 outside, 0 onbound, +1 inside
    statv = view(stat, ivec)

    tol = max(abs(rtol * lbar), abs(atol))
    inpoly2!(vert, node, edge, flip, tol, statv)
    stat
end

"""
    inpoly2_mat(vert, node, edge, fTol, stats)

INPOLY2_MAT the local m-code version of the crossing-number
test. Loop over edges; do a binary-search for the first ve-
rtex that intersects with the edge y-range; do crossing-nu-
mber comparisons; break when the local y-range is exceeded.
"""
function inpoly2!(vert, node, edge, flip::Bool, veps::AbstractFloat, stat::AbstractVector{<:Integer})
    nvrt = size(vert, 1) # the points to be checked
    nnod = size(node, 1) # the vertices of the polygon
    nedg = size(edge, 1) # the indices of the vertices connected
    
    ix = flip + 1
    iy = 2 - flip
    #----------------------------------- loop over polygon edges
    for epos = 1:nedg

        inod = edge[epos,1]  # from
        jnod = edge[epos,2]  # to
        if node[inod,iy] > node[jnod,iy]
            inod, jnod = jnod, inod
        end

        #------------------------------- calc. edge bounding-box
        xone = node[inod,ix]
        yone = node[inod,iy]
        xtwo = node[jnod,ix]
        ytwo = node[jnod,iy]

        xmin = min(xone, xtwo) - veps
        xmax = max(xone, xtwo) + veps
        ymin = yone - veps
        ymax = ytwo + veps

        ydel = ytwo - yone
        xdel = xtwo - xone
        feps = veps * hypot(xdel, ydel)

        #------------------------------- find top vert[:,iy] < ymin
        ilow = 1
        iupp = nvrt
        while ilow < iupp - 1    # binary search
            imid = ilow + (iupp-ilow) ÷ 2
            if vert[imid,iy] < ymin
                ilow = imid
            else
                iupp = imid
            end
        end
        while ilow > 0 && vert[ilow,iy] >= ymin
            ilow = ilow - 1
        end

        #------------------------------- calc. edge-intersection
        # loop over all points with y ∈ [ymin,ymax]
        for jpos = ilow+1:nvrt
            stat[jpos] == 0 && continue
            ypos = vert[jpos,iy]
            ypos > ymax && break 
            xpos = vert[jpos,ix]

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

#=
# adaptation to `inpolygon` interface
struct Engwirda <: MembershipCheckAlgorithm end
export Engwirda

function inpoly(p::AbstractArray{T}, poly, ::Union{Engwirda,Type{Engwirda}}) where T <:Real
    vert = [[p[1] p[2]];]
    inpoly2a(vert, poly)
end
function inpolygon(pp, poly, ::Union{Engwirda,Type{Engwirda}})
    xpos = [p[1] for p in pp]
    ypos = [p[2] for p in pp]
    inpoly2a([xpos ypos], poly)
end
function inpoly2a(vert, poly)
    node, edge = convertpoly(poly)
    stat, bnds = inpoly2(vert, node, edge)
    convertres(stat, bnds)
end

function convertpoly(poly::AbstractVector)
    x = [p[1] for p in poly]
    y = [p[2] for p in poly]
    convertpoly([x y])
end
function convertpoly(node::AbstractMatrix)
    n = size(node, 1)
    edge = [[1:n...] [2:n...;1]]
    node, edge
end

function convertres(stat::AbstractVector{Bool}, bnds::AbstractVector{Bool})
    n = length(stat)
    [ b ? 0 : s ? 1 : -1 for (s, b) in zip(stat, bnds)]
end

function inpolygon1(pp::AbstractVector, poly, algo = HaoSun())
    inpolygon.(pp, Ref(poly), Ref(algo))
end
function inpolygon1(pp::AbstractMatrix{T}, poly, algo = HaoSun()) where T<:Real
    if size(pp, 1) == 1
        inpolygon([pp[1,1]; pp[1,2]], poly, algo)
    else
        p = [[pp[k,1]; pp[k,2]] for k in axes(pp,1)]
        inpolygon.(p, Ref(poly), Ref(algo))
    end
end
=#
