module FittingFunction


export PL
export SBPL



"""
    PL(E,N,β;E₀ = 1.)

Computes a power-law with spectral index 'β' and normalization 'N' at input 'E'. You can normalize the power-law at 'E₀'.


# Examples
```jldoctest
PL(3.,1.,-1.)

# output

0.3333333333333333
```
"""
function PL(E,N,β;E₀ = 1.)
    spec = .^(E./E₀,β)
    return N .* spec
end





"""
    SBPL(E,N,α₁,α₂,bE)

Computes a smoothly joint broken power-law with spectral index '$\alpha_1$' and '$\alpha_2$', a break at 'Eb' and normalization 'N' at input 'E'.


# Examples
```jldoctest
SBPL(6.,1.,-1.,1.5,5.)

# output

0.26290682760247974
```
""
function SBPL(E,N,α₁,α₂,Eb)
    f1 = PL(E,N,α₁)
    n = PL(Eb,N,α₁) ./ PL(Eb,1.0,α₂)
    f2 = PL(E,n,α₂)
    return ifelse.(E .<= Eb, f1, f2)
end



end
