
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
function polydemo(id::Integer=1; args...)
    [demo1, demo2, demo3, demo4][id](;args...)
end

#-----------------------------------------------------------
function demo1(;tol=1e-2, r=10^5, do_plot=true, which=3)
    do_plot && println("""
        INPOLY2 provides fast point-in-polygon queries for ob-\n 
        jects in R^2. Like INPOLYGON, it detects points inside\n
        and on the boundary of polygonal geometries.\n
""")

    node = [ 2.0 0; 8 4; 4 8; 0 4; 5 0; 7 4; 3 6; 6 1 # inner nodes
             3 3; 5 3; 5 5; 3 5 ]                     # outer nodes  
    edge = [ 1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 1
             9 10; 10 11; 11 12; 12 9 ]

    dx = 10 * sqrt(1/r)
    a = [(x,y) for x in -1:dx:9 for y in -1:dx:9]
    xpos, ypos = [p[1] for p in a], [p[2] for p in a]

    stat, bnds = inpoly2([xpos ypos], node, edge, atol=tol)
    if do_plot
        p = plot(title="demo1 (r $r, atol $tol)", legend=:topleft)
        p = plotpoints(p, xpos, ypos, stat, bnds, which)
        p = plotpolygon(p, node, edge)
        display(p)
    end
    stat, bnds
end

#-----------------------------------------------------------
function demo2(;r=2500, tol=1e-3, do_plot=true, which=3)
    do_plot && println("""
        INPOLY2 supports multiply-connected geometries, consi-\n
        sting of arbitrarily nested sets of outer + inner bou-\n
        ndaries.\n
""")

    xpos, ypos, node, edge = testdata("lakes.msh", r)
    stat, bnds = inpoly2([xpos ypos], node, edge, rtol=tol) 
    if do_plot
        p = plot(title="demo2 - lakes (r $r, rtol $tol)", legend=:topleft)
        p = plotpolygon(p, node, edge)
        p = plotpoints(p, xpos, ypos, stat, bnds, which)
        display(p)
    end
    stat, bnds
end

#-----------------------------------------------------------
function demo3(;r=2500, tol=1e-3, do_plot=true, which=3)
    do_plot && println("""
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
    stat = @time inpolygon1(xvert, xnode)
    =#
    println("inpoly2")
    stat, bnds = @time inpoly2([xpos ypos], node, edge, rtol=tol)
    if do_plot
        p = plot(title="demo3 - coast (r $r, rtol $tol)", legend=:topleft)
        p = plotpolygon(p, node, edge)
        p = plotpoints(p, xpos, ypos, stat, bnds, which)
        display(p)
    end
    stat, bnds
end

#-----------------------------------------------------------
function demo4(;r=5*10^5, tol=1e-2, do_plot=true, which=3)
    do_plot && println("""
        INPOLY2 provides fast point-in-polygon queries for ob-\n 
        jects in R^2. Like INPOLYGON, it detects points inside\n
        and on the boundary of polygonal geometries.\n
    """)

    node = [0.0 0; 3 4; 7 5; 10 0]
    edge = [1 2; 2 3; 3 4; 4 1]

    rpts = testbox(r, [-1.0 -1.0], [ 11.0 6.0])
    xpos, ypos = rpts[:,1], rpts[:,2]

    stat, bnds = inpoly2(rpts, node, edge, atol=tol)
    if do_plot
        p = plot(title="demo4 (r $r, atol $tol)", legend=:topleft)
        p = plotpoints(p, xpos, ypos, stat, bnds, which)
        p = plotpolygon(p, node, edge)
        display(p)
    end
    stat, bnds
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
        p = plot!(p, xpos, ypos, color=:purple, width=1.5)
    end
    p
end

function plotpoints(p, xpos, ypos, inside, bnds, which::Integer=3)
    ms = 1.0
    outside = (!).(inside)
    onbound = bnds .& ((inside .& (which & 1 == 1)) .| (outside .& (which & 2 == 2)))
    scatter!(p, xpos[onbound], ypos[onbound], markershape=:x, markersize=ms, color=:black)
    scatter!(p, xpos[inside], ypos[inside], markershape=:cross, markersize=ms, color=:blue)
    scatter!(p, xpos[outside], ypos[outside], markershape=:cross, markersize=ms, color=:red)
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

