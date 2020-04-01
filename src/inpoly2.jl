"""
    function inpoly2(vert, node, edge[, fTOL])

Check all points defined by `vert` are inside, outside or on bounds of
polygon defined by `node` and `edge`.
`vert` and `node` are matrices with the x-coordinates in the first
and the y-coordinates in the second column.
`edge` is an integer matrix of same size a `node`. The first column contains
the index of the starting node, the second column the index of the ending node.
The polygon needs to be closed. It may contain unconnected cycles.
"""
function inpoly2(vert::AbstractMatrix{T}, node::AbstractMatrix{T}, edge::AbstractMatrix{<:Integer}=Int[], fTOL::T=eps(T)^0.85) where T<:Real

    nnod = size(node,1)
    nvrt = size(vert,1)
    if isempty(edge)
        edge = [(1:nnod-1)  (2:nnod); nnod 1]
    else
        edge = copy(edge)
    end
    node = copy(node)

    if size(node,2) != 2 || size(edge,2) != 2 || size(vert,2) != 2
        throw(ArgumentError("inpoly2:incorrectDimensions"))
    end

    if minimum(edge) < 1 || maximum(edge) > nnod
        throw(ArgumentError("inpoly2:invalidInputs: Invalid EDGE input array."))
    end

#-------------- flip to ensure the y-axis is the "long" axis
    vmin = minimum(vert,dims=1)
    vmax = maximum(vert,dims=1)
    ddxy = vmax - vmin

    lbar = sum(ddxy) / 2

    if ddxy[1] > ddxy[2]
        vert[:,:] = vert[:,[2,1]]
        node[:,:] = node[:,[2,1]]
    end

#----------------------------------- sort points via y-value
    swap = node[edge[:,2],2] .< node[edge[:,1],2]
    edge[swap,:] = edge[swap,[2,1]]

    ivec = sortperm(view(vert,:,2))
    vert = view(vert,ivec,:)

    stat, bnds = inpoly2_mat(vert, node, edge, fTOL, lbar)
    stat[ivec] = stat
    bnds[ivec] = bnds
    stat, bnds
end

"""
    inpoly2_mat(vert, node, edge, fTol, lbar)

INPOLY2_MAT the local m-code version of the crossing-number
test. Loop over edges; do a binary-search for the first ve-
rtex that intersects with the edge y-range; do crossing-nu-
mber comparisons; break when the local y-range is exceeded.
"""
function inpoly2_mat(vert, node, edge, fTOL, lbar)

    veps = fTOL * lbar

    nvrt = size(vert, 1) # the points to be checked
    nnod = size(node, 1) # the vertices of the polygon
    nedg = size(edge, 1) # the indices of the vertices connected

    stat = falses(nvrt)
    bnds = falses(nvrt)

    #----------------------------------- loop over polygon edges
    for epos = 1:nedg

        inod = edge[epos,1]  # from
        jnod = edge[epos,2]  # to

        #------------------------------- calc. edge bounding-box
        yone = node[inod,2]
        ytwo = node[jnod,2]
        xone = node[inod,1]
        xtwo = node[jnod,1]

        xmin = min(xone, xtwo) - veps
        xmax = max(xone, xtwo) + veps

        ymin = yone - veps # assumption yone <= ytwo
        ymax = ytwo + veps

        ydel = ytwo - yone
        xdel = xtwo - xone
        feps = veps * hypot(xdel, ydel)

        #------------------------------- find top vert[:,2] < ymin
        ilow = +1
        iupp = nvrt
        while ilow < iupp - 1    # binary search
            imid = ilow + (iupp-ilow) ÷ 2
            if vert[imid,2] < ymin
                ilow = imid
            else
                iupp = imid
            end
        end
        if vert[ilow,2] >= ymin
            ilow = ilow - 1
        end

        #------------------------------- calc. edge-intersection
        for jpos = ilow+1:nvrt
            bnds[jpos] && continue

            xpos = vert[jpos,1]
            ypos = vert[jpos,2]

            if ypos <= ymax
                if xpos >= xmin
                    if xpos <= xmax
                    #------------------- compute crossing number
                    mul1 = ydel * (xpos - xone)
                    mul2 = xdel * (ypos - yone)

                        if feps >= abs(mul2 - mul1)
                            #------------------- BNDS -- approx. on edge
                            bnds[jpos]= true
                            stat[jpos]= true

                        elseif mul1 < mul2
                            if ypos >= yone && ypos < ytwo
                                #------------------- left of edge
                                stat[jpos] = ~stat[jpos]
                            end
                        end
                    end
                else
                    if (ypos >= yone && ypos <  ytwo)
                        #------------------- advance crossing number
                        stat[jpos] = ~stat[jpos]
                    end
                end
            else
                break # done -- due to the sort
            end
        end
    end
    stat, bnds
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
