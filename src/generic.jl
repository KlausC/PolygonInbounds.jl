"""
    PolygonType(nodes[, edges])

Polygon given by `nodes` and `edges`. Several formats are supported:
 * nodes as nx2 matrix
 * nodes as n-vector of 2-vectors or 2-tuples of x-y-coordinates
 * edges as nx2 matrix of indices
 * edges as n-vector of indices`
 * edges empty

Incomplete edges are filles up by best guess. The columns must be permutations of `1:n`.
"""
struct PolygonType{U,E}
    nodes::U
    edges::E
    PolygonType(n::U, ed::E) where {U,E} = new{U,E}(n, ed)
    PolygonType(n::U) where U = new{U,Nothing}(n, nothing)
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
    see `InOnOut` and `InOnBool`
"""
abstract type AbstractOutputFormat end
"""
    InOnOut{IN,ON,OUT}

Expected format of output is one of three constant values.
`InOnOut{1,0,-1}` returns `1`, `0`, `-1` if point is in, on, or out respectively.
"""
struct InOnOut{IN,ON,OUT} <: AbstractOutputFormat end

"""
    InOnBool

Expected output is a pair of booleans, the first one, `inout`, indicating the point
is inside, the second one, `bnds`, indicating the point is close to the border.
Depending on the used algorithm, the value of `inout` may be invalid, if `bnds` is true.
"""
struct InOnBool <: AbstractOutputFormat end

# interface for obove structures standard cases
function Base.length(poly::PolygonType)
    size(poly.nodes, 1)
end
function edge(poly::PolygonType{<:Any,Union{Nothing,<:AbstractArray{<:Any,0}}}, i::Integer)
    n = length(poly)
    j = i ÷ n + 1
    i, j
end
function edge(poly::PolygonType{<:Any,<:AbstractVector}, i::Integer, n::Integer)
    n == 1 ? i : poly.edges[i]
end
function edge(poly::PolygonType{<:Any,<:AbstractMatrix}, i::Integer, n::Integer)
    poly.edges[i,n]
end

function vertex(poly::PolygonType{<:AbstractMatrix{<:Real}}, v::Integer, xy::Integer)
    poly.nodes[v,xy]
end
function vertex(poly::PolygonType{<:AbstractVector}, v::Integer, xy::Integer)
    poly.nodes[v][xy]
end

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
    sortperm(p.points[:,ix])
end

function convertout(to::Type{T}, from::Type{S}, args...) where {S<:AbstractOutputFormat,T<:AbstractOutputFormat}
    T == S && return length(args) == 1 ? args[1] : args
    convertto(convertfrom(from, args...))
end

function convertfrom(::Type{InOnBool}, stat::T, bnds::T) where T<:Union{BitVector,Bool}
    stat, bnds
end
function convertfrom(::Type{InOnBool}, stat::T, bnds::T) where T<:AbstractVector{Bool}
    convert(BitVector, stat), convert(BitVector, bnds)
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
function convertfrom(::Type{OT}, stat::AbstractVector{<:Integer}) where OT<:InOnOut
    n = length(stat)
    inout = BitArray(undef, n)
    bnds = BitArray(undef, n)
    for i in 1:n
        a, b = convertfrom(OT, stat[i])
        inout[i] = a
        bnds[i] = b
    end
    inout, bnds
end

function convertto(::Type{<:InOnBool}, stat::T, bnds::T) where T<:Union{BitVector,Bool}
    stat, bnds
end
function convertto(::Type{InOnOut{IN,ON,OUT}}, stat::Bool, bnds::Bool) where {IN,ON,OUT}
    bnds ? ON : stat ? IN : OUT
end
function convertto(::Type{OT}, stat::T, bnds::T) where {OT<:InOnOut,T<:AbstractVector{Bool}}
    convertto.(OT, stat, bnds) 
end

function checkedge(edge)
    if N == 0
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
    nothing
end
