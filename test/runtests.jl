using FittingFunction
using Test

@testset "FittingFunction.jl" begin
    # Write your tests here.
    # PL
    @test PL(3.,1.,-1.) == 0.3333333333333333
    #
    # SBPL
    @test SBPL(3.,1.,-1.,-1.5,5.) == 0.3333333333333333
    @test SBPL(6.,1.,-1.,-1.5,5.) == 0.15214515486254615
    #
    @test Band(300.,-2.,-1.,200.,1.) == -0.45304697140984085
    #
    @test XAbs([0.3,0.5]) == [0.7322879401681527, 0.9290803992951834]
    #
    @test Extinction(5500.,1;gal="MW",Rv=FFGals["MW"],z=0.) == 0.06017150537641227
end
