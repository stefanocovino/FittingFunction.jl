module FittingFunction

using DustExtinction



export Band
export Extinction
export FFGals
export PL
export SBPL
export SBPL2
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

Computes a power-law with spectral index 'α' and normalization 'N' at input 'E'. You can normalize the power-law at 'E₀'.


# Examples
```jldoctest
CPL(3.,1.,-1.,2.)

# output

0.6693904804452895
```
"""
function CPL(E,N,α,Ec; E0=1.)
    return N .* .^(E./E0,α) .* exp.(-E ./ Ec)
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
    PL(E,N,α;E0=1.)

Computes a power-law with spectral index 'α' and normalization 'N' at input 'E'. You can normalize the power-law at 'E₀'.


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
    return map(e -> XAbsorption(e, NH=NH, z=z), E)
end


end
