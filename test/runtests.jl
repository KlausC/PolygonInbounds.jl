using Test
using PolygonInbounds
using .Demo

@testset "PolygonInboundsTests" begin
@testset "argument checks" begin
    
    @test inpoly2([0.0 0], [0.0 0]) == [false true]
    @test_throws ArgumentError inpoly2([0.0 0], [0.0 0], [1;2])
    @test inpoly2([0.0 0], [-1.0 1; 1 -1]) == [false true]
    @test inpoly2([0.0 0], [-1.0 1; 1 -1], [1; 2]) == [false false]
    @test inpoly2([0.0 0], [-1.0 1; 1 -1], [2; 1]) == [false true]
    @test_throws ArgumentError inpoly2([0.0 0], [-1.0 1; 1 -1], [2 1])
    @test inpoly2([0.0 0], [-1.0 1; 1 -1], [1 1]) == [false false]
    @test_throws ArgumentError inpoly2([0.0 0], [-1.0 1; 1 -1], [2 1; 2 2])
    @test_throws ArgumentError inpoly2([0.0 0], [-1.0 1; 1 -1], [2 2; 1 2])
    @test_throws ArgumentError  inpoly2([0.0 0], [-1.0 1; 1 -1], [1 1 1 2 3 4])
end

let stat = polydemo(4, r=10^5, tol=1e-2, do_plot=false)
    @testset "functionality" begin
        @test 300 <= count(stat[:,2]) < 600
        @test 29000 <= count(stat[:,1]) < 32000
        @test 67000 <= count((!).(stat[:,1])) < 70000
    end
end

@testset "utilities" begin
    @test PolygonInbounds.mean_overlap(0.0, 0.0, 2.0, 8.0, 10.0) == 0
    @test PolygonInbounds.mean_overlap(2.0, 0.0, 2.0, 8.0, 10.0) == 1.5
    @test PolygonInbounds.mean_overlap(8.0, 0.0, 2.0, 8.0, 10.0) == 6
    @test PolygonInbounds.mean_overlap(11.0, 0.0, 2.0, 8.0, 10.0) == 6
end

@testset "means estimations" begin
    r = 100
    nodes = [0.0 0; 3 4; 7 5; 10 0]
    edges = [1 2; 2 3; 3 4; 4 1]
    poly = PolygonInbounds.PolygonMesh(nodes, edges)
    rpts = PolygonInbounds.PointsInbound(Demo.testbox(r, [-1.0 -1.0], [ 11.0 6.0]))
    x = PolygonInbounds.estimate_effort(rpts, [-1.0 -1.0], [11.0 6.0], poly, [0.0 0.0], [10.0 5.0], 1e-6)
    @test length(x) == 4
end

end

# to improve test coverage
for i = 1:5
    polydemo(i)
end

