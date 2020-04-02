using Test
using PolygonInbounds

@testset "argument checks" begin
    
    @test inpoly2([0.0 0], [0.0 0]) == [0]
    @test_throws ArgumentError inpoly2([0.0 0], [0.0 0], [1;2])
    @test inpoly2([0.0 0], [-1.0 1; 1 -1]) == [0]
    @test inpoly2([0.0 0], [-1.0 1; 1 -1], [1; 2]) == [-1]
    @test inpoly2([0.0 0], [-1.0 1; 1 -1], [2; 1]) == [0]
    @test_throws ArgumentError inpoly2([0.0 0], [-1.0 1; 1 -1], [2 1])
    @test_throws ArgumentError inpoly2([0.0 0], [-1.0 1; 1 -1], [1 1])
    @test_throws ArgumentError inpoly2([0.0 0], [-1.0 1; 1 -1], [2 1; 2 2])
    @test_throws ArgumentError inpoly2([0.0 0], [-1.0 1; 1 -1], [2 2; 1 2])
end

let stat = polydemo(4, r=10^5, tol=1e-2, do_plot=false)
    @testset "functionality" begin
        @test 300 <= count(stat .== 0) < 600
        @test 29000 <= count(stat .== 1) < 32000
        @test 67000 <= count(stat .== -1) < 70000
    end
end
