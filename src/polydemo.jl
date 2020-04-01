
using Plots
export polydemo

"""
%POLYDEMO run a series of point-in-polygon demos.
%   POLYDEMO(II) runs the II-th demo, where +1 <= II <= +3. 
%   Demo problems illustrate varying functionality for the 
%   INPOLY2 routine.
%
%   See also INPOLY2, INPOLYGON

%   Darren Engwirda : 2018 --
%   Email           : darren.engwirda@columbia.edu
%   Last updated    : 28/10/2018
"""
function polydemo(id::Integer=1)
    [demo1, demo2, demo3, demo4][id]()
end

function demo1()
    println("""
        INPOLY2 provides fast point-in-polygon queries for ob-\n 
        jects in R^2. Like INPOLYGON, it detects points inside\n
        and on the boundary of polygonal geometries.\n
""")

    node = [
        2.0 0       # outer nodes
        8 4
        4 8
        0 4
        5 0
        7 4
        3 6
        6 1
        3 3         # inner nodes
        5 3
        5 5
        3 5
        ] ;
    edge = [
        1 2         # outer edges
        2 3
        3 4
        4 5
        5 6
        6 7
        7 8
        8 1
        9 10        # inner edges
        10 11
        11 12
        12 9
        ] ;

    a = [(x,y) for x in -1:0.2:9 for y in -1:0.2:9]
    xpos, ypos = [p[1] for p in a], [p[2] for p in a]

    stat, bnds = inpoly2([xpos ypos], node, edge)
    p = plotpolygon(plot(title="demo1"), node, edge)
    p = plotpoints(p, xpos, ypos, stat, bnds)
    display(p)
end

function demo2(r=2500)
#-----------------------------------------------------------
    println("""
        INPOLY2 supports multiply-connected geometries, consi-\n
        sting of arbitrarily nested sets of outer + inner bou-\n
        ndaries.\n
""")

    xpos, ypos, node, edge = testdata("lakes.msh", r)
    stat, bnds = inpoly2([xpos ypos], node, edge) 
    p = plotpolygon(plot(title="demo2 - lakes"), node, edge)
    p = plotpoints(p, xpos, ypos, stat, bnds)
    display(p)
end

function demo3(r=2500)
    println("""
        INPOLY2 implements a "pre-sorted" variant of the cros-\n
        sing-number test - returning queries in approximately \n
        O((N+M)*LOG(N)) time for configurations consisting of \n
        N test points and M edges. This is often considerably \n
        faster than conventional approaches that typically re-\n
        quire O(N*M) operations.\n
        """);

    xpos, ypos, node, edge = testdata("coast.msh", r)
    
    #=
    println("inpolynom")
    xvert = [zip(xpos, ypos)...]
    xnode = [(node[k,:] for k = axes(node, 1))...]
    stat, bnds = @time inpolygon1(xvert, xnode)
    =#
    println("inpoly2")
    stat, bnds = @time inpoly2([xpos ypos], node, edge)

    p = plotpolygon(plot(title="demo3 - coast"), node, edge)
    p = plotpoints(p, xpos, ypos, stat, bnds)
    display(p)
end

function demo4()
    println("""
        INPOLY2 provides fast point-in-polygon queries for ob-\n 
        jects in R^2. Like INPOLYGON, it detects points inside\n
        and on the boundary of polygonal geometries.\n
    """)

    node = [0.0 0; 3 4; 7 4; 7 0]
    edge = [1 2; 2 3; 3 4; 4 1]

    rpts = testbox(10^5, [-1.0 -1.0], [ 8.0 5.0])
    xpos, ypos = rpts[:,1], rpts[:,2]

    stat, bnds = inpoly2(rpts, node, edge, 0.05)
    p = plot(title="demo4")
    p = plotpoints(p, xpos, ypos, stat, bnds)
    p = plotpolygon(p, node, edge)
    display(p)
end

#----- Utility functions

function testdata(dataset, r::Integer)
    filepath = @__DIR__
    node, edge = loadmsh(joinpath(filepath, "..", "test-data", dataset))

    emid = .5 * node[edge[:,1],:] + .5 * node[edge[:,2],:]

    ma, mi =  maximum(node, dims=1), minimum(node, dims=1)
    rpts = testbox(r, mi, ma)

    vert = r <= 100 ? [rpts;] : [node; emid; rpts]
    vert[:,1], vert[:,2], node, edge
end

function testbox(r::Integer, mi, ma)
    mi = reshape(mi, 1, 2)
    ma = reshape(ma, 1, 2)
    half = (ma + mi) / 2
    scal = ma - mi
    rpts = (rand(r, 2) .- 0.5).* scal .* 1.1 .+ half
end

function plotpolygon(p, node::Matrix{T}, edge::Matrix{<:Integer}) where T<:Real
    d = Dict((edge[k,1], edge[k,2]) for k in 1:size(edge,1))
    while !isempty(d)
        start = j = minimum(keys(d))
        n = length(d)
        xpos = Vector{T}(undef, n+1)
        ypos = Vector{T}(undef, n+1)
        k = 0
        while k <= n
            k += 1
            xpos[k] = node[j,1]
            ypos[k] = node[j,2]
            j == start && k > 1 && break
            j = pop!(d, j)
        end
        resize!(xpos, k)
        resize!(ypos, k)
        p = plot!(p, xpos, ypos, color=:purple, width=2)
    end
    p
end

function plotpoints(p, xpos, ypos, stat, bnds)
    ms = 1.0
    scatter!(p, xpos[map(!,stat)], ypos[map(!,stat)], markershape=:cross, markersize=ms, color=:red)
    scatter!(p, xpos[stat], ypos[stat], markershape=:cross, markersize=ms, color=:blue)
    scatter!(p, xpos[bnds], ypos[bnds], markershape=:x, markersize=ms, color=:black)
    display(p)
    p
end

function loadmsh(file)
    nodes = Float64[]
    edges = Int[]
    m = 0
    n = 0
    open(file) do io
        name =""
        while !eof(io)
            line = readline(io)
            startswith(line, "#") && continue
            ls = split(line, "=")
            if length(ls) == 2
                name = ls[1]
                if name == "point"
                    n = parse(Int, ls[2])
                    nodes = Vector{Float64}(undef, 2n)
                    m = 0
                elseif name == "edge2"
                    n = parse(Int, ls[2])
                    edges = Vector{Int}(undef, 2n)
                    m = 0
                end
                continue
            end
            ls = split(line, ";")
            if length(ls) >= 2
                m >= 2n - 1 && continue
                if name == "point"
                    nodes[m+=1] = parse(Float64, ls[1])
                    nodes[m+=1] = parse(Float64, ls[2])
                elseif name == "edge2"
                    edges[m+=1] = parse(Int, ls[1]) + 1
                    edges[m+=1] = parse(Int, ls[2]) + 1
                end
            end
        end
    end
    res(nodes) = permutedims(reshape(nodes, 2, length(nodes)÷2))
    res(nodes), res(edges)
end

