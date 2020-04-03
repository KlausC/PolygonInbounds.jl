"""
    PolygonMesh(nodes[, edges])

Polygon given by `nodes` and `edges`. Several formats are supported:
 * nodes as nx2 matrix
 * nodes as n-vector of 2-vectors or 2-tuples of x-y-coordinates
 * edges as nx2 matrix of indices
 * edges as n-vector of indices`
 * edges empty

Incomplete edges are filles up by best guess. The columns must be permutations of `1:n`.
"""
struct PolygonMesh{A,U,E<:Union{Nothing,AbstractArray{<:Integer}}}
    nodes::U
    edges::E
    function PolygonMesh(n::U, ed::E) where {U,E}
        A = ed isa AbstractMatrix ? max(size(ed, 2) - 2, 0) : 0
        # the number of stored area indices per edge
        checkarguments(n, ed)
        new{A,U,E}(n, ed)
    end
    PolygonMesh(n::U) where U = PolygonMesh(n, nothing)
end

"""
    PointsInbound(p)

Vector of points. Several formats are supported for points:
 * as 2-vector of real numbers (set of single point)
 * as nx2 matrix of real numbers - each row represents `(x, y)`
 * as n-vector of 2-vectors or 2-tuples of x-y-coordinates
"""
struct PointsInbound{V}
    points::V
    PointsInbound(p::V) where V = new{V}(p)
end

"""
    see `InOnOut` and `InOnBit`
"""
abstract type AbstractOutputFormat end
"""
    InOnOut{IN,ON,OUT}

Expected format of output is one of three constant values.
`InOnOut{1,0,-1}` returns `1`, `0`, `-1` if point is in, on, or out respectively.
"""
struct InOnOut{IN,ON,OUT} <: AbstractOutputFormat end

"""
    InOnBit

Expected output is a pair of booleans, the first one, `inout`, indicating the point
is inside, the second one, `bnds`, indicating the point is close to the border.
Vectors are stored in `BitArray`s of 3 dimensions.

Depending on the used algorithm, the value of `inout` may be invalid, if `bnds` is true.
"""
struct InOnBit <: AbstractOutputFormat end

# interface for above structures standard cases

# access functions for PolygonMesh

function nodecount(poly::PolygonMesh)
    size(poly.nodes, 1)
end
function edgecount(poly::PolygonMesh)
    size(poly.edges, 1)
end

function edgeindex(poly::PolygonMesh{0,<:Any,<:Union{Nothing,AbstractArray{<:Any,0}}}, i::Integer, n::Integer)
    n == 1 ? i : i % nodecount(poly) + 1
end
function edgeindex(poly::PolygonMesh{<:Any,<:Any,<:AbstractVector}, i::Integer, n::Integer)
    n == 1 ? i : poly.edges[i]
end
function edgeindex(poly::PolygonMesh{<:Any,<:Any,<:AbstractMatrix}, i::Integer, n::Integer)
    poly.edges[i,n]
end

function areaindex(poly::PolygonMesh{A}, i::Integer, a::Integer) where A
    A >= 1 ? poly.edges[i,a+2] : 1
end
function areacount(poly::PolygonMesh{A}) where A
    if A >= 1
        mi, ma = extrema(view(poly.edge, :, 3:A))
        mi >= 0 || throw(ArgumentError("negative area index detected"))
        ma
    else
        1
    end
end

function vertex(poly::PolygonMesh{<:Any,<:AbstractMatrix{<:Real}}, v::Integer, xy::Integer)
    poly.nodes[v,xy]
end
function vertex(poly::PolygonMesh{<:Any,<:AbstractVector}, v::Integer, xy::Integer)
    poly.nodes[v][xy]
end

# access functions for PointsInbounds

function Base.length(p::PointsInbound{<:AbstractVector{<:Real}})
    1
end
function Base.length(p::PointsInbound{<:AbstractMatrix{<:Real}})
    size(p.points, 1)
end
function Base.length(p::PointsInbound{<:AbstractVector})
    length(p.points)
end

function vertex(p::PointsInbound{<:AbstractVector{<:Real}}, v::Integer, xy::Integer)
    p.points[xy]
end
function vertex(p::PointsInbound{<:AbstractMatrix{<:Real}}, v::Integer, xy::Integer)
    p.points[v,xy]
end
function vertex(p::PointsInbound{<:AbstractVector}, v::Integer, xy::Integer)
    p.points[v][xy]
end

# standard functions for PointsInbound

function Base.minimum(p::PointsInbound)
    n = length(p)
    minimum([[vertex(p, k, 1) for k = 1:n] [vertex(p, k, 2) for k = 1:n]], dims = 1)
end
function Base.maximum(p::PointsInbound)
    n = length(p)
    maximum([[vertex(p, k, 1) for k = 1:n] [vertex(p, k, 2) for k = 1:n]], dims = 1)
end
function Base.sortperm(p::PointsInbound, ix::Integer)
    n = length(p)
    sortperm([vertex(p, k, ix) for k in 1:n])
end
function Base.sortperm(p::PointsInbound{<:AbstractMatrix}, ix::Integer)
    sortperm(view(p.points, :, ix))
end


# conversion between different output formats

function convertout(to::Type{T}, from::Type{S}, args...) where {S<:AbstractOutputFormat,T<:AbstractOutputFormat}
    T == S && return length(args) == 1 ? args[1] : args
    convertto(to, convertfrom(from, args...)...)
end

function convertfrom(::Type{<:InOnBit}, stat::T, bnds::T) where T<:Bool
    stat, bnds
end
function convertfrom(::Type{<:InOnBit}, stat::T) where T<:AbstractArray{Bool}
    stat
end
function convertfrom(::Type{InOnOut{IN,ON,OUT}}, stat::Integer) where {IN,ON,OUT}
    if stat == IN
        true, false
    elseif stat == OUT
        false, false
    elseif stat == ON
        true, true
    else
        false, true
    end
end
function convertfrom(::Type{OT}, stat::Union{AbstractVector{T},AbstractMatrix{T}}) where {T<:Integer,OT<:InOnOut}
    n, m = size(stat)
    inout = BitArray(undef, n, 2, m)
    for j = 1:m
        for i in 1:n
            a, b = convertfrom(OT, stat[i,j])
            inout[i,1,j] = a
            inout[i,2,j] = a
        end
    end
    inout, bnds
end

function convertto(::Type{<:InOnBit}, stat::T, bnds::T) where T<:Bool
    stat, bnds
end
function convertto(::Type{<:InOnBit}, stat::T, bnds::T) where T<:BitArray
    stat
end
function convertto(::Type{InOnOut{IN,ON,OUT}}, stat::Bool, bnds::Bool) where {IN,ON,OUT}
    bnds ? ON : stat ? IN : OUT
end
function convertto(::Type{OT}, stat::T) where {OT<:InOnOut,T<:AbstractArray{Bool}}
    convertto.(OT, view(stat, :,1,:), view(stat, :, 2, :))
end

# check functions for PolygonMesh constructor

function checknodes(nodes::AbstractMatrix)
    size(nodes, 2) >= 2 || throw(ArgumentError("nodes require at least 2 coordinates"))
    nothing
end
function checknodes(nodes::AbstractVector)
    nothing
end
function checknodes(nodes)
    nothing
end

function checkarguments(nodes, edges::Union{Nothing,AbstractArray{<:Integer,0}})
    checknodes(nodes)
    nothing
end
function checkarguments(nodes, edges::AbstractVector{<:Integer})
    checknodes(nodes)
    n = size(nodes, 1)
    mi, ma = extrema(edges)
    0 < mi <= ma <= n || throw(ArgumentError("invalid edge indices"))
    nothing
end
function checkarguments(nodes, edges::AbstractMatrix{<:Integer})
    checknodes(nodes)
    n, m = size(nodes)
    m >= 2 || throw(ArgumentError("edges matrix requires at least 2 columns"))
    m <= 5 || throw(ArgumentError("edges matrix accepts at most 3 area columns"))
    mi, ma = extrema(view(edges, :, 1:2))
    0 < mi <= ma <= n || throw(ArgumentError("invalid edge indices"))
    if m > 2
        mi, ma = extrema(view(edges, :, 3:m))
        0 <= mi <= 4096 || throw(ArgumentError("invalid area indices ($mi, $ma)"))
    end
    if m <= 2 || m == 3 && mi == ma
        check(i) = isperm(view(edges, :, i)) ||
        throw(ArgumentError("edges column $i is no permutation"))
        check(1)
        check(2)
    end
    nothing
end
