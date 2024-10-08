using DataFrames
using FittingFunction
using Test


@testset "FittingFunction.jl" begin
    # Write your tests here.
    @test PL(3.,1.,-1.) == 0.3333333333333333
    #
    @test SBPL(3.,1.,-1.,-1.5,5.) == 0.3333333333333333
    @test SBPL(6.,1.,-1.,-1.5,5.) == 0.15214515486254615
    #
    @test Band(300.,-2.,-1.,200.,1.) == -0.45304697140984085
    #
    @test XAbs([0.3,0.5]) == [0.7322879401681527, 0.9290803992951834]
    #
    @test Extinction(5500.,1;gal="MW",Rv=FFGals["MW"],z=0.) == 0.06017150537641227
    #
    @test SBPL2([4.,6.,8.],1.,-1.,-1.5,-2.,5.,7.) == [0.25, 0.06804138174397717, 0.04133986423538423]
    #
    @test CPL(3.,1.,-1.,2.) == 0.0743767200494766
    #
    @test Gaussian(3.,1.,0.1,2.) == 4.009419869036418
    #
    @test GaussAbs(2.,1.,0.1,2.) == 0.018510395162776326
    #
    @test Counts2Mag(100.,10.) == (20.0, 0.10857362047581294)
    #
    @test Mag2Counts(20.,0.1) == (100.0, 9.210340371976184)
    #
    @test Stokes2Pol(1.,0.1,0.1,0.) == (1.0, 0.14142135623730953, 0.39269908169872414, 0.0, 0.0, 0.0, 0.0)
    #
    @test Pol2Stokes(1.,0.1,0.1,0.) == (1.0, 0.09800665778412417, 0.019866933079506124, 0.0, 0.0, 0.0, 0.0)
    #
    @test FittingFunction.VoigtFunctTG(1e-5,1.) == 0.36788094737611143
    #
    @test TauVoigt(1494.5 .* 1e-8,1e18,1e5,0.,[1495.05*1e-8,0.54,8.106e8]) == 0.5408331361420329
    #
    @test typeof(GetAtomicData()) == typeof(DataFrame())
    #
    x = [4.,6.,8.,1.,3.,5.,20.]
    mask = SigmaClip(x)
    @test x[mask] == [4.,6.,8.,1.,3.,5.]
    #
    @test Jy2Ecm2sA([1e-3,2e-3,3e-3],[5000,5500,6000]) == [1.2000000000000002e-15,1.9834710743801656e-15,2.5e-15]
    #
    @test Ecm2sA2Jy([1.2e-15,2e-15,2.5e-15],[5000,5500,6000]) == [0.001,0.002016666666666667,0.0029999999999999996]
    #
    @test FourierPeriodogram([1.,2.,3.,4.],1.) == ([0.0, 0.25], [100.0, 8.000000000000002])
    #
end
