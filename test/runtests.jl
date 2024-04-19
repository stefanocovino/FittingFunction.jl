using FittingFunction
using Test

@testset "FittingFunction.jl" begin
    # Write your tests here.
    # PL
    @test PL(3.,1.,-1.) == 0.3333333333333333
    #
    @test SBPL(3.,1.,-1.) == 0.3333333333333333
    @test SBPL(6.,1.,-1.) == 0.26290682760247974
end
