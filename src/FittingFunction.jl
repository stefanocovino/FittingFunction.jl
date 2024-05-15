module FittingFunction

using CSV
using DataFrames
using DustExtinction
using PhysicalConstants.CODATA2018
using Unitful




export Band
export CPL
export Counts2Mag
export Extinction
export FFGals
export GaussAbs
export Gaussian
export GetAtomicData
export Mag2Counts
export PL
export Pol2Stokes
export SBPL
export SBPL2
export Stokes2Pol
export TauVoigt
export XAbs



"""
    Band(E,α,β,E₀,A)

Compute the so-called 'Band function', at input energy 'E'. 'α' is the low-energy index, 'β' the high-energy index and 'E₀' the break energy. 'A' is the normalization. The peak-energy is 'Ep = (2+α)E₀'.

# Examples


```jldoctest
Band(300.,-2.,-1.,200.,1.)

# output

-0.45304697140984085
```
"""
function Band(E,α,β,E₀,A)
    f1 = A .* (E ./ 100.).^α .* exp.(-E ./ E₀)
    f2 = A .* (E ./ 100.).^β .* exp.(β .- α) .* ((α .- β) .* E₀ ./ 100.).^(α .- β)
    return ifelse.(E .<= (α .- β) .* E₀, f1, f2)
end



"""
    CPL(E,N,α,Ec;E0=1.)

Computes a cut-off power-law with spectral index 'α' and normalization 'N' at input 'E' and cut-off energy 'Ec'. You can normalize the power-law at 'E0'.


# Examples
```jldoctest
CPL(3.,1.,-1.,2.)

# output

0.0743767200494766
```
"""
function CPL(E,N,α,Ec; E0=1.)
    return N .* .^(E./E0,α) .* exp.(-E ./ Ec)
end



"""
    Counts2Mag(cts,ects;zp=25.0,ezp=0.0)

Convert counts (or flux) 'cts', with uncertainty 'ects', to magnitude scale with zero-point 'zp' and relative uncertainty.


# Examples
```jldoctest
Counts2Mag(100,10)

# output

(20.0, 0.10857362047581294)
```
"""
function Counts2Mag(cts,ects;zp=25.0,ezp=0.)
    mag = -2.5 .* log10.(cts) .+ zp
    emag = (2.5 ./ log(10)) .* (ects ./ cts)
    emag = sqrt.(emag .^2 .+ ezp .^ 2)
    return mag, emag
end





FFGals = Dict("MW" => 3.12, "SMC" => 2.74, "LMC" => 3.41, "Any" => 4.0)
"""
    Extinction(wave,EBV;gal="SMC",Rv=FFGals["SMC"],z=0.)

Compute the UV/optical/extinction at the input wavelength 'wave', in Angstrom. 'EBV' is E(B-V) in magnitudes, 'gal' is one of the galaxy extinction recipes listed in the 'FFGals' (exported) dictionary. 'Rv' is the seelctive extinction and 'z' is the redshift of the absorpber. References about the adopted extinction curves are discussed in the documentation of the [DustExtinction](https://juliaastro.org/DustExtinction.jl/stable/) package.


# Examples


```jldoctest
Extinction(5500.,1,gal="MW",Rv=FFGals["MW"],z=0.)

# output

0.06017150537641227
```
"""
function Extinction(wave,EBV;gal="SMC",Rv=FFGals["SMC"],z=0.)
    # rest frame
    if gal == "SMC"
        adlr = G03_SMCBar().(wave ./ (1 .+ z)) .* Rv * EBV
    elseif gal == "LMC"
        adlr = G03_LMCAve().(wave ./ (1 .+ z)) .* Rv * EBV
    elseif gal == "MW"
        adlr = G16(Rv=Rv,f_A=1).(wave ./ (1 .+ z)) .* Rv * EBV
    else
        adlr = G16(Rv=Rv,f_A=1).(wave ./ (1 .+ z)) .* Rv * EBV
    end
    absr = .^(10, -adlr ./ 2.5)
    #
    return absr
end






"""
    GaussAbs(E,Ed,σ,El)

Computes a Gaussian absorption with depth 'Ed' at energy 'E' with center in 'El' and width (FWHM/2.35) 'σ'. 'A.


# Examples
```jldoctest
GaussAbs(2.,1.,0.1,2.)

# output

0.018510395162776326
```
"""
function GaussAbs(E,Ed,σ,El)
    return exp.(-(Ed / sqrt(2*pi)) .* (1 ./ σ) .* exp.((E .- El).^2 ./ 2*σ.^2))
end







"""
    Gaussian(E,A,σ,El)

Computes a Gaussian at energy 'E' with center in 'El' and width (FWHM/2.35) 'σ'. 'A is the normalization.


# Examples
```jldoctest
Gaussian(3.,1.,0.1,2.)

# output

4.009419869036418
```
"""
function Gaussian(E,A,σ,El)
    return A .* (1/sqrt(2*pi)) * (1 ./ σ) .* exp.((E .- El).^2 ./ 2*σ.^2)
end





AvailableAtomicTables = ["FitLyman","JitrikBunge04","VALD3"]

"""
    GetAtomicData(table::String="")::DataFrame

Provides a table with atomic data among those available. Calling the function with no parameter shows the available tables. At present we have:

* "FitLyman". This is the original table generated for the ["FitLyman" ESO-MIDAS](https://ui.adsabs.harvard.edu/abs/1995Msngr..80...37F/abstract) package.
* "JitrikBunge04". This is a collection of atomic data prepared by [Jitrik & Bunge (2004)](https://ui.adsabs.harvard.edu/abs/2004JPCRD..33.1059J/abstract).
* "VALD03". This is a collection of atomic data extracted from the [VALD3](http://vald.astro.uu.se/) database.



# Examples
```julia
fly = GetAtomicData("FitLyman")
```
"""
function GetAtomicData(table::String="")::DataFrame
    if table == ""
      for t in AvailableAtomicTables
        println("Available atomic data tables:")
        println("\t"*t)
      end
      return DataFrame()
    elseif table in AvailableAtomicTables
      tb = DataFrame(CSV.File(joinpath(pkgdir(FittingFunction),"data",table*".csv")))
      return tb
    else
      println("Warning! Table not recognized.")
      return DataFrame()
    end
end






"""
    Mag2Counts(mag,emag;zp=25.0)

Convert magnitudes 'mag' with uncertainty 'emag' to counts (or flux) with zero-point 'zp'.


# Examples
```jldoctest
Mag2Counts(20,0.1)

# output

(100.0, 9.210340371976184)
```
"""
function Mag2Counts(mag, emag; zp=25.0)
    cts = 10 .^ (-0.4 .* (mag .- zp))
    ects = (emag .* cts) ./ (2.5 ./ log(10.0))
    return cts, ects
end



"""
    PL(E,N,α;E0=1.)

Computes a power-law with spectral index 'α' and normalization 'N' at input 'E'. You can normalize the power-law at 'E0'.


# Examples
```jldoctest
PL(3.,1.,-1.)

# output

0.3333333333333333
```
"""
function PL(E,N,α; E0=1.)
    spec = .^(E./E0,α)
    return N .* spec
end






"""
    Pol2Stokes(i, p, t, c; ep=0.0, et=0.0, ec=0.0)

Computes the Stokes parameters for a given polarisation degree and position angle. 'i' is total intensity, 'p' the polarisation degree, 't' the position angle (randians) and 'c' the circular polarisation. The errors 'ep', 'ew, and 'e' are optional.
The ouput is: intensity 'q', 'u', 'v' Stokes parameters and respective errors.


# Examples
```jldoctest
Pol2Stokes(1.,0.1,0.1,0.)

# output

(1.0, 0.09800665778412417, 0.019866933079506124, 0.0, 0.0, 0.0, 0.0)
```
"""
function Pol2Stokes(i, p, t, c; ep=0.0, et=0.0, ec=0.0)
    q = p.*i.*cos.(2 .*t).*cos.(2 .*c)
    u = p.*i.*sin.(2 .*t).*cos.(2 .*c)
    v = p.*i.*sin.(2 .*c)
    #
    eq = sqrt.((q.*ep./p).^2 .+ (2 .*u.*et).^2 .+ (2 .*v.*cos.(2 .*t).*ec).^2)
    eu = sqrt.((u.*ep./p).^2 .+ (2 .*q.*et).^2 .+ (2 .*v.*sin.(2 .*t).*ec).^2)
    ev = sqrt.((v.*ep./p).^2 .+ (2 .*i.*p.*ec.*cos.(2 .*c)).^2)
    return i, q, u, v, eq, eu, ev
end






"""
    SBPL(E,N,α,β,Eb)

Computes a smoothly joint broken power-law with spectral index 'α' and 'β', a break at 'Eb', normalization 'N' and input energy'E'.


# Examples
```jldoctest
SBPL(6.,1.,-1.,-1.5,5.)

# output

0.15214515486254615
```
"""
function SBPL(E,N,α,β,Eb)
    f1 = PL(E,N,α)
    n = PL(Eb,N,α) ./ PL(Eb,1.0,β)
    f2 = PL(E,n,β)
    return ifelse.(E .<= Eb, f1, f2)
end



"""
    SBPL2(E,N,α,β,γ,Eb1,Eb2)

Computes a smoothly joint double broken power-law with spectral index 'α', 'β' and 'γ', and two breaks at 'Eb1' and 'Eb2'. The normalization is 'N' and the input energy 'E'.


# Examples
```jldoctest
SBPL2([4.,6.,8.],1.,-1.,-1.5,-2.,5.,7.)

# output

3-element Vector{Float64}:
 0.25
 0.06804138174397717
 0.04133986423538423
```
"""
function SBPL2(E,N,α,β,γ,Eb1,Eb2)
    f1 = PL(E,N,α)
    n12 = PL(Eb1,N,α) ./ PL(Eb1,1.0,α)
    f2 = PL(E,n12,β)
    n23 = PL(Eb2,n12,β) ./ PL(Eb2,1.0,γ)
    f3 = PL(E,n23,γ)
    return ifelse.(E .<= Eb1, f1, ifelse.(E .<= Eb2, f2, f3))
end



"""
    Stokes2Pol(i, q, u, v; eq=0.0, eu=0.0, ev=0.0)

Computes the polarisation degree and position angle given the Stokes parameters: 'i', total intensity, and 'q', 'u' and 'v'. The errors 'eq', 'eu, and 'ev' are optional.
The ouput is: intensity 'i', polarisation degree, 'p', position angle (randians) 'theta' and circular polarisation 'chi', with respective errors.


# Examples
```jldoctest
Stokes2Pol(1.,0.1,0.1,0.)

# output

(1.0, 0.14142135623730953, 0.39269908169872414, 0.0, 0.0, 0.0, 0.0)
```
"""
function Stokes2Pol(i, q, u, v; eq=0.0, eu=0.0, ev=0.0)
    pl = sqrt.(q.^2 .+ u.^2) ./ i
    pol = sqrt.(q.^2 .+ u.^2 .+ v.^2) ./ i
    theta = 0.5 .* atan.(u,q)
    chi = 0.5 .* atan.(v, pl)
    #
    epl = sqrt.((q.*eq).^2 .+ (u.*eu).^2)./pl
    epol = sqrt.((q.*eq).^2 .+ (u.*eu).^2 .+ (v.*ev).^2)./pol
    etheta = 0.5.*sqrt.((u.*eq).^2 .+ (q.*eu).^2)./(pl.^2)
    echi = 0.5.*sqrt.((v.*epl).^2 .+ (pl.*ev).^2)./(pol.^2)
    return i, pol, theta, chi, epol, etheta, echi
end




c0 = ustrip(uconvert(u"cm"*u"s^-1",SpeedOfLightInVacuum))
em = ustrip(uconvert(u"g",ElectronMass))
ee = ustrip(ElementaryCharge)/3.3356e−10;


"""
    TauVoigt(λ,NHI,DopplerBroadening,z,transition)

Computes the opacity due to a given absorption transition. ``\\lambda`` are the wavelenghths (``cm``), ``N`` the column density (``cm^{-2}``), 'DopplerBroadening' is expressed in (``cm~s^{-1}``), 'z' is the source redshift and 'transition' is a tuple formed by the central wawelength (``cm``), the oscillator strength and the damping coefficient (``s^{-1}``) for the given transition. The Doppler broadening factor is typicaly due to either an intrinsic thermal broadening (``b_K = \\sqrt{2KT/m}``) or to a turbulent motion of the gas (``b_T = \\sqrt{\\sigma_T}``, where ``\\sigma_T`` is the inner velocity dispersion). These two factors can be summed in quadrature, i.e. ``b = \\sqrt{b_K^2 + b_T^2}``.


# Examples
```jldoctest
TauVoigt(range(start=1494.5,stop=1495.5,step=0.1) .* 1e-8,1e18,1e5,0.,[1495.05*1e-8,0.54,8.106e8])

# output

11-element Vector{Float64}:
  0.5408331361420329
  0.8079604230648749
  1.3357686994197793
  2.618871809580839
  7.282358478732964
 66.40911467182524
 66.40911467182524
  7.282358478732964
  2.618871809580839
  1.3357686994184856
  0.8079604230648749
```
"""
function TauVoigt(λ,N,DopplerBroadening,z,transition)
    LambdCentr, OscillStrength, DampCoeff = transition
    #
    coeffNum = 4 * sqrt(pi^3) * ee^2
    coeffDem = em * c0
    a = (LambdCentr * DampCoeff) / (4 * pi * DopplerBroadening)
    fact = a * (coeffNum/coeffDem) * (OscillStrength/DampCoeff)
    v = c0 * (LambdCentr .- λ ./ (1 .+ z)) / (LambdCentr * DopplerBroadening)
    #
    res = VoigtFunctTG(a,v)
    #
    return N .* res .* fact
end




"""
    VoigtFuncTG(a,v)

Approximates the Voigt (or line broadening function) function. It was introduced by [Tepper-Garcia (2006)](https://ui.adsabs.harvard.edu/abs/2006MNRAS.369.2025T/abstract). It is of [high accuracy](https://en.wikipedia.org/wiki/Voigt_profile) provided that ``a \\le 10^{-4}``. The ``a`` and ``v`` parameters are defined as:
``a = \\frac{\\Gamma \\lambda_c}{4\\pi b 10^{13}}`` and ``v = \\frac{(\\lambda_c - \\lambda)c}{b \\lambda_c \\sqrt{2 \\log 2}}``, where ``\\lambda_c`` is the central wavelength of the transition, ``c`` the speed of light, ``b`` the Doppler broadening factor and ``\\Gamma`` the damping coefficient.


# Examples
```jldoctest
FittingFunction.VoigtFunctTG(1e-5,1.)

# output

0.36788094737611143
```
"""
function VoigtFunctTG(a,v)
    F = v.^2
    H0 = exp.(-F)
    Q = 1.5 ./ F
    res1 = H0
    res2 = a ./ (pi.^0.5 .* F)
    res3 = H0.^2 .* (4 .* F.^2 .+ 7 .* F .+ 4 .+ Q) .- Q .- 1
    return res1 .- res2 .* res3
end







function XAbsorption(E; NH=1e20, z=0)
    ren = E .* (1 + z)
    if ren < 0.03
        return 0.0
    elseif ren > 10.0
        return 1.0
    end
    if NH < 0
        return 0.0
    end
    if 0.03 <= ren < 0.1
        C1 = 17.3
        C2 = 608.1
        C3 = -2150.0
    elseif 0.1 <= ren < 0.284
        C1 = 34.6
        C2 = 267.9
        C3 = -476.1
    elseif 0.284 <= ren < 0.400
        C1 = 78.1
        C2 = 18.8
        C3 = 4.3
    elseif 0.4 <= ren < 0.532
        C1 = 71.4
        C2 = 66.8
        C3 = -51.4
    elseif 0.532 <= ren < 0.707
        C1 = 95.5
        C2 = 145.8
        C3 = -61.1
    elseif 0.707 <= ren < 0.867
        C1 = 308.9
        C2 = -380.6
        C3 = 294.0
    elseif 0.867 <= ren < 1.303
        C1 = 120.6
        C2 = 169.3
        C3 = -47.7
    elseif 1.303 <= ren < 1.840
        C1 = 141.3
        C2 = 146.8
        C3 = -31.5
    elseif 1.840 <= ren < 2.471
        C1 = 202.7
        C2 = 104.7
        C3 = -17.0
    elseif 2.471 <= ren < 3.210
        C1 = 342.7
        C2 = 18.7
        C3 = 0.0
    elseif 3.210 <= ren < 4.038
        C1 = 352.2
        C2 = 18.7
        C3 = 0.0
    elseif 4.038 <= ren < 7.111
        C1 = 433.9
        C2 = -2.4
        C3 = 0.75
    elseif 7.111 <= ren < 8.331
        C1 = 629.0
        C2 = 30.9
        C3 = 0.0
    elseif ren >= 8.331
        C1 = 701.2
        C2 = 25.2
        C3 = 0.0
    end
    xcs = 1e-24 .* (C1 .+ C2 .* ren .+ C3 .* ren.^2) ./ ren.^3
    return exp.(-(xcs .* NH))
end



"""
    XAbs(E; NH=1e20, z=0)

Computes the effective absorption cross section per hydrogen atom following [Morrison & McCammon (1983)](https://ui.adsabs.harvard.edu/abs/1983ApJ...270..119M/abstract). 'E' is the energy (KeV, in the 0.03-10 KeV range), 'NH' is the column density of hydrogen atom (particle per square cm) and 'z' is the redshift.


# Examples
```jldoctest
XAbs([0.5,1.25,2.], NH=1e20, z=0)

# output

3-element Vector{Float64}:
 0.9290803992951834
 0.9868927382232576
 0.9957079871273042
```
"""
function XAbs(E; NH=1e20, z=0)
    return map(e -> FittingFunction.XAbsorption(e, NH=NH, z=z), E)
end


end
